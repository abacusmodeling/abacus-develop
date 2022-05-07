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
printf "%20s %7s %8s %8s %6s %6s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n" \
"example" "Natoms" "EneCut" "k-points" "NProc" "Niter" "TotTime" "1stSCF" "SCF/iter" \
"s_c%" "c_bands%" "s_bands%" "h_psi%" "vloc%" "vnl%" "FFT%" "stress%" "force%"> $outf

for i in `cat $allcase`;do
    if [[ ! -f $i/result.log ]];then
        echo "Error: NO FILE $i/result.log"
        continue
        #printf "%20s" $i >> sum.dat 
    fi
    
    basis=`awk '$1=="basis_type"{print $2}' ${i}/INPUT | tr [A-Z] [a-z]`
    #echo $basis
    if [[ "$basis" == "pw" ]];then
        natoms=`sed -n '/ELEMENT NATOM/,/----/'p ${i}/result.log| sed  '1d;$d' | awk 'BEGIN{a=0}{a+=$2}END{print a}'`
        encut=`grep "energy cutoff for wavefunc" ${i}/OUT.*/running*|awk '{print $NF}'`
        kpt=`grep -A 1 "KPOINTS" ${i}/result.log | tail -1 | awk '{print $2}'`
        nproc=`grep -A 1 "KPOINTS" ${i}/result.log | tail -1 | awk '{print $3}'`
        niter=`sed -n '/ITER   ETOT(eV)/,/><><><><>/'p ${i}/result.log | wc -l|awk '{print $1-2}'`
        tottime=`awk '$1=="total"{printf"%.2f", $2}' ${i}/result.log`
        scf1=`grep -A 1 "ITER   ETOT(eV)" ${i}/result.log | awk 'END{printf"%.2f", $NF}'`
        totalscf=`awk '$2=="self_consistent"{print $3}' ${i}/result.log`
        scfpiter=`awk -v a=$totalscf -v b=$scf1 -v c=$niter 'BEGIN{printf"%.2f",(a-b)/(c-1)}'`
        fft=`awk '$2=="FFT3D"{printf"%.1f",$6}' ${i}/result.log`
        hpsi=`awk '$2=="h_psi"{printf"%.1f",$6}' ${i}/result.log`
        vloc=`awk '$2=="vloc"{printf"%.1f",$6}' ${i}/result.log`
        vnl=`awk '$2=="vnl"{printf"%.1f",$6}' ${i}/result.log`
        sc=`awk '$2=="self_consistent"{printf"%.1f",$6}' ${i}/result.log`
        cbands=`awk '$2=="c_bands"{printf"%.1f",$6}' ${i}/result.log`
        sbands=`awk '$2=="sum_band"{printf"%.1f",$6}' ${i}/result.log`
        stress=`awk '$2=="cal_stress"{printf"%.1f",$6}' ${i}/result.log`
        force=`awk '$2=="cal_force_nl"{printf"%.1f",$6}' ${i}/result.log`
    elif [[ "$basis" == "lcao" ]];then
        natoms=`sed -n '/ELEMENT ORBITALS/,/----/'p ${i}/result.log| sed  '1d;$d' | awk 'BEGIN{a=0}{a+=$4}END{print a}'`
        encut=`grep "energy cutoff for wavefunc" ${i}/OUT.*/running*|awk '{print $NF}'`
        kpt=`grep -A 1 "KPOINTS" ${i}/result.log | tail -1 | awk '{print $2}'`
        nproc=`grep -A 1 "KPOINTS" ${i}/result.log | tail -1 | awk '{print $3}'`
        niter=`sed -n '/ITER   ETOT(eV)/,/><><><><>/'p ${i}/result.log | wc -l|awk '{print $1-2}'`
        tottime=`awk '$1=="total"{printf"%.2f", $2}' ${i}/result.log`
        scf1=`grep -A 1 "ITER   ETOT(eV)" ${i}/result.log | awk 'END{printf"%.2f", $NF}'`
        totalscf=`awk '$1=="ELEC_scf"{print $3}' ${i}/result.log`
        scfpiter=`awk -v a=$totalscf -v b=$scf1 -v c=$niter 'BEGIN{printf"%.2f",(a-b)/(c-1)}'`
        fft="-"
        hpsi="-"
        vloc=`awk '$2=="vlocal"{printf"%.1f",$6}' ${i}/result.log`
        vnl="-"
        sc=`awk '$1=="ELEC_scf"{printf"%.1f",$6}' ${i}/result.log`
        cbands=`awk '$2=="cal_bands"{printf"%.1f",$6}' ${i}/result.log`
        sbands=`awk '$2=="sum_bands"{printf"%.1f",$6}' ${i}/result.log`
        stress=`awk '$2=="evaluate_vl_stress"{printf"%.1f",$6}' ${i}/result.log`
        force=`awk '$2=="evaluate_vl_force"{printf"%.1f",$6}' ${i}/result.log`
    else
        echo "ERROR: UNKNOW basis type $basis"
        continue
    fi

    printf "%20s %7s %8s %8s %6s %6s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s %8s\n" \
    $i $natoms $encut $kpt $nproc $niter $tottime $scf1 $scfpiter $sc $cbands $sbands $hpsi $vloc $vnl $fft \
    $stress $force >> $outf

done
