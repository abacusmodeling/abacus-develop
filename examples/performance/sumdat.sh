#!/bin/bash

allcase=allcase
outf=sum.dat

if [[ $# > 0 ]];then
    allcase=$1
    if [[ $# > 1 ]];then
        outf=$2
    fi
fi

test -f $outf && rm $outf

#title
printf "%20s %15s %7s %8s %8s %6s %6s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n" \
"example" "ks_solver" "Natoms" "EneCut" "k-points" "NProc" "Niter" "TotTime" "1stSCF" "SCF/iter" \
"Run%" "h_psi%" "vloc%" "vnl%" "stress%" "force%" "MaxResSize" > $outf

for i in `cat $allcase`;do
    if [[ ! -f $i/result.log ]];then
        echo "Error: NO FILE $i/result.log"
        continue
        #printf "%20s" $i >> sum.dat 
    fi
    
    basis=`awk '$1=="basis_type"{print $2}' ${i}/INPUT | tr [A-Z] [a-z]`
    solver=`grep ks_solver ${i}/OUT.*/INPUT | awk '{print $2}'`
    #echo $basis
    if [[ "$basis" == "pw" ]];then
        natoms=`sed -n '/ELEMENT NATOM/,/----/'p ${i}/result.log| sed  '1d;$d' | awk 'BEGIN{a=0}{a+=$2}END{print a}'`
        encut=`grep "energy cutoff for wavefunc" ${i}/OUT.*/running*|awk '{print $NF}'`
        kpt=`grep -A 1 "KPOINTS" ${i}/result.log | tail -1 | awk '{print $2}'`
        nproc=`grep -A 1 "KPOINTS" ${i}/result.log | tail -1 | awk '{print $3}'`
        niter=`grep "ELEC=" $i/OUT.*/running* | awk -F "=" '{print $3}'| awk -F "-" 'END{print $1}'`
        tottime=`awk '$1=="total"{printf"%.2f", $2}' ${i}/result.log`
        scf1=`awk '{if($1 ~ /^CG1$|^DAV1$|^GE1$|^GV1$/) {printf"%.2f", $NF}}' ${i}/result.log`
        totalscf=`awk '$2=="Run"{print $3}' ${i}/result.log`
        scfpiter=`awk 'BEGIN{printf"%.2f",(ARGV[1]-ARGV[2])/(ARGV[3]-1)}' $totalscf $scf1 $niter`
        #fft=`awk '$2=="FFT3D"{printf"%.1f",$6}' ${i}/result.log`
        hpsi=`awk '$2=="h_psi"{printf"%.1f",$6}' ${i}/result.log`
        vloc=`awk '$2=="vloc"{printf"%.1f",$6}' ${i}/result.log`
        vnl=`awk '$2=="vnl"{printf"%.1f",$6}' ${i}/result.log`
        sc=`awk '$2=="Run"{printf"%.1f",$6}' ${i}/result.log`
        #cbands=`awk '$2=="c_bands"{printf"%.1f",$6}' ${i}/result.log`
        #sbands=`awk '$2=="sum_band"{printf"%.1f",$6}' ${i}/result.log`
        stress=`awk '$2=="cal_stress"{printf"%.1f",$6}' ${i}/result.log`
        force=`awk '$2=="cal_force_nl"{printf"%.1f",$6}' ${i}/result.log`
    elif [[ "$basis" == "lcao" ]];then
        natoms=`sed -n '/ELEMENT ORBITALS/,/----/'p ${i}/result.log| sed  '1d;$d' | awk 'BEGIN{a=0}{a+=$4}END{print a}'`
        encut=`grep "energy cutoff for wavefunc" ${i}/OUT.*/running*|awk '{print $NF}'`
        kpt=`grep -A 1 "KPOINTS" ${i}/result.log | tail -1 | awk '{print $2}'`
        nproc=`grep -A 1 "KPOINTS" ${i}/result.log | tail -1 | awk '{print $3}'`
        niter=`grep "ELEC=" $i/OUT.*/running* | awk -F "=" '{print $3}'| awk -F "-" 'END{print $1}'`
        tottime=`awk '$1=="total"{printf"%.2f", $2}' ${i}/result.log`
        scf1=`awk '{if($1 ~ /^CG1$|^DAV1$|^GE1$|^GV1$/) {printf"%.2f", $NF}}' ${i}/result.log`
        totalscf=`awk '$2=="Run"{print $3}' ${i}/result.log`
        scfpiter=`awk 'BEGIN{printf"%.2f",(ARGV[1]-ARGV[2])/(ARGV[3]-1)}' $totalscf $scf1 $niter`
        #fft="-"
        hpsi="-"
        vloc=`awk '$2=="vlocal"{printf"%.1f",$6}' ${i}/result.log`
        vnl="-"
        sc=`awk '$2=="Run"{printf"%.1f",$6}' ${i}/result.log`
        #cbands=`awk '$2=="cal_bands"{printf"%.1f",$6}' ${i}/result.log`
        #sbands=`awk '$2=="sum_bands"{printf"%.1f",$6}' ${i}/result.log`
        stress=`awk '$2=="evaluate_vl_stress"{printf"%.1f",$6}' ${i}/result.log`
        force=`awk '$2=="evaluate_vl_force"{printf"%.1f",$6}' ${i}/result.log`
    else
        echo "ERROR: UNKNOW basis type $basis"
        continue
    fi
    maxres=`grep "Maximum resident set size" ${i}/time.log | awk '{print $NF}'`

    if [[ $solver == "" ]]; then solver="-";fi
    if [[ $natoms == "" ]]; then natoms="-";fi
    if [[ $encut == "" ]]; then encut="-";fi
    if [[ $kpt == "" ]]; then kpt="-";fi
    if [[ $nproc == "" ]]; then nproc="-";fi
    if [[ $niter == "" ]]; then niter="-";fi
    if [[ $tottime == "" ]]; then tottime="-";fi
    if [[ $scf1 == "" ]]; then scf1="-";fi
    if [[ $scfpiter == "" ]]; then scfpiter="-";fi
    if [[ $sc == "" ]]; then sc="-";fi
    if [[ $hpsi == "" ]]; then hpsi="-";fi
    if [[ $vloc == "" ]]; then vloc="-";fi
    if [[ $vnl == "" ]]; then vnl="-";fi
    if [[ $stress == "" ]]; then stress="-";fi
    if [[ $force == "" ]]; then force="-";fi
    if [[ $maxres == "" ]]; then maxres="-";fi

    printf "%20s %15s %7s %8s %8s %6s %6s %8s %8s %8s %8s %8s %8s %8s %8s %8s %s\n" \
    $i $solver $natoms $encut $kpt $nproc $niter $tottime $scf1 $scfpiter $sc $hpsi $vloc $vnl \
    $stress $force $maxres >> $outf

done
