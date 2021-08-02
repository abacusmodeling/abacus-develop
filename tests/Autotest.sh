#!/bin/bash

# ABACUS executable path
abacus=abacus
# number of mpi processes
np=1
# threshold with unit: eV
threshold=0.0000001
# check accuracy
ca=8
# regex of case name
case="^[^#].*01_PW_.*$"
# enable AddressSanitizer
sanitize=false

while getopts a:n:t:c:s:r,g flag
do
    case "${flag}" in
        a) abacus=${OPTARG};;
        n) np=${OPTARG};;
		t) threshold=${OPTARG};;
		c) ca=${OPTARG};;
		s) sanitize=${OPTARG};;
		r) case=${OPTARG};;
		g) g=true;; #generate test reference
    esac
done

echo "-----AUTO TESTS OF ABACUS ------"
echo "ABACUS path: $abacus";
echo "Number of cores: $np";
echo "Test accuracy: $threshold eV"
echo "Check accuaracy: $ca"
echo "Test cases: $case"
echo "Generate reference: $g"
echo "--------------------------------"
echo ""

#----------------------------------------------------------
# define a function named 'check_out'
#----------------------------------------------------------
check_out(){
	#------------------------------------------------------
	# input file $1 is 'result.out' in each test directory
	#------------------------------------------------------
	outfile=$1

	#------------------------------------------------------
	# outfile = result.out
	#------------------------------------------------------
	properties=`awk '{print $1}' $outfile`

	#------------------------------------------------------
	# jd = job description
	#------------------------------------------------------
	if test -e "jd"; then
		jd=`cat jd`
 		echo " [  ------  ] $jd"
	fi

	#------------------------------------------------------
	# check every 'key' word
	#------------------------------------------------------
	for key in $properties; do

		#--------------------------------------------------
		# calculated value
		#--------------------------------------------------
		cal=`grep "$key" result.out | awk '{printf "%.'$ca'f\n",$2}'`

		#--------------------------------------------------
		# reference value
		#--------------------------------------------------
		ref=`grep "$key" result.ref | awk '{printf "%.'$ca'f\n",$2}'`

		#--------------------------------------------------
		# computed the deviation between the calculated
		# and reference value
		#--------------------------------------------------
		deviation=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$ca'f\n",x-y}'`
		deviation1=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$ca'f\n",y-x}'`

		#--------------------------------------------------
		# computed the deviation between the calculated
		# and reference value for descriptors in DeePKS
		#--------------------------------------------------
		if [ $key == "descriptor" ]; then
			check_file descriptor.dat
			state=`echo $?`
			if [ $state == "0" ]; then
				let failed++
				break
			fi
		fi

		if [ $key == "jle" ]; then
			check_file jle.orb
			state=`echo $?`
			if [ $state == "0" ]; then
				let failed++
				break
			fi
		fi

		if [ $key == "totaltimeref" ]; then
			# echo "time=$cal ref=$ref"
			break
		fi		


		#--------------------------------------------------
		# If deviation < threshold, then the test passes,
		# otherwise, the test prints out warning
		# Daye Zheng found bug on 2021-06-20,
		# deviation should be positively defined
		#--------------------------------------------------
		if [ ! -n "$deviation" ]; then
            echo -e "\e[1;31m Fatal Error! :(  \e[0m"
			let failed++
			break
        else
			if [ $(echo "sqrt($deviation*$deviation) < $threshold"|bc) = 0 ]; then
				echo -e "\e[1;31m [  FAILED  ] \e[0m"\
					"$key cal=$cal ref=$ref deviation=$deviation"
				let failed++
				break
			else
				#echo "$key cal=$cal ref=$ref deviation=$deviation"
				#echo "[ PASS ] $key"
				echo -e "\e[1;32m [      OK  ] \e[0m $key"
			fi
		fi
		let ok++
	done
}

#------------------------------------------------------------
# define a function named 'check_file' to compare every data 
#------------------------------------------------------------
check_file(){
	#--------------------------------------------------
	# input file $1 is 'descriptor.dat' or 'jle.orb'
	#--------------------------------------------------
	if [ $1 == "jle.orb" ]; then
		outfile=OUT.autotest/$1
	else
		outfile=$1
	fi
	reffile=$1.ref

	#---------------------------------------------
	# compare every data
	#---------------------------------------------
	row_ref=$(wc -l < $reffile)
	row_cal=$(wc -l < $outfile)
	if [ $row_ref -ne $row_cal ]; then
		echo -e "\e[1;31m [  FAILED  ] \e[0m"\
			"$1  the number of data rows is different !"
		return 0
	fi

	for irow in `seq 1 $row_ref`; do

		if [ $(( $irow % 100 )) == "0" ]; then
			echo "$1  Compare No.$irow row ~"
		fi

		sum1=`awk -F " " 'NR=="'"$irow"'" {print NF}' $reffile`
		sum2=`awk -F " " 'NR=="'"$irow"'" {print NF}' $outfile`

		if [ $sum1 -ne $sum2 ]; then
			echo -e "\e[1;31m [  FAILED  ] \e[0m"\
				"$1  the number of datas in row $irow is different !"
			return 0
		fi

		if [ $sum1 == "0" ]; then
			continue
		fi

		for num in `seq 1 $sum1`; do
			ref=`awk -F " " 'NR=="'"$irow"'" {print $"'"$num"'"}' $reffile`
			cal=`awk -F " " 'NR=="'"$irow"'" {print $"'"$num"'"}' $outfile`
			if [ $ref != $cal ]; then
				deviation=`awk 'BEGIN {x='$ref';y='$cal';printf "%.'$ca'f\n",x-y}'`
				if [ $(echo "sqrt($deviation*$deviation) < $threshold"|bc) = 0 ]; then
					echo -e "\e[1;31m [  FAILED  ] \e[0m"\
						"$1  row=$irow column=$num cal=$cal ref=$ref deviation=$deviation"
					return 0
				fi
    		fi
		done
	done

	return 1
}


#---------------------------------------------
# the file name that contains all of the tests
#---------------------------------------------

test -e CASES || echo "Plese specify tests"
test -e CASES || exit 0
which $abacus > /dev/null || echo "Error! ABACUS path was wrong!!"
which $abacus > /dev/null || exit 0

testdir=`cat CASES | grep -E $case`
failed=0
ok=0
report=""
repo="$(realpath ..)/"

if [ "$sanitize" == true ]; then
	echo "Testing with Address Sanitizer..."
	mkdir ../html
	echo -e "# Address Sanitizer Diagnostics\n" > ../html/README.md
	report=$(realpath ../html/README.md)
fi

for dir in $testdir; do
	cd $dir
	echo -e "\e[1;32m [  RUN     ]\e[0m $dir"
	TIMEFORMAT=' [  ------  ] Time elapsed: %R seconds'
	#parallel test
	time {
		if [ "$sanitize" == true ]; then
			ASAN_OPTIONS="log_path=asan" mpirun -np $np $abacus > log.txt
			echo -e "## Test case ${dir}\n" >> ${report}
			for diagnostic in asan.*; do
				echo -e "### On process id ${diagnostic}\n" >> ${report}
				echo -e "\`\`\`bash" >> ${report}
				cat ${diagnostic} >> ${report}
				echo -e "\`\`\`\n" >> ${report}
			done
		else
			mpirun -np $np $abacus > log.txt
		fi
		#$abacus > log.txt
		test -d OUT.autotest || echo "Some errors happened in ABACUS!"
		test -d OUT.autotest || exit 0
		if test -z $g
		then
			../tools/catch_properties.sh result.out
			check_out result.out
		else
			../tools/catch_properties.sh result.ref
		fi
	}
	echo ""
	cd ../
done

if [ "$sanitize" == true ]; then
	if [[ `uname` == "Darwin" ]]; then
		sed -i '' "s,${repo},,g" ${report}
	else
		sed -i "s,${repo},,g" ${report}
	fi
fi

if [ -z $g ]
then
if [ $failed -eq 0 ]
then
	echo -e "\e[1;32m [  PASSED  ] \e[0m $ok test cases passed."
else
	echo -e "\e[1;31m [  FAILED  ] \e[0m $failed test cases out of $[ $failed + $ok ] failed."
	exit 1
fi
else
echo "Generate test cases complete."
fi
