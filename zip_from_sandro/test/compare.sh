#!/bin/bash

file=$1
#/Users/sandro/Programs/Phaistos/phaistos/build_debug/test/test_ff ${file}
#dir=~/Programs/Phaistos/phaistos/ff_scripts/test
rm \#*
rm ${dir}/md.log
rm gromacs_ene.dat 
rm charmm_ene.dat


pdb2gmx -f ${file} -o ${dir}/temp.gro <<EOF 
8
6
EOF
exit
editconf -bt cubic -f ${dir}/temp.gro -o ${dir}/temp_box.gro -c -d 1.0  
grompp -f ${dir}/mini.mdp -c ${dir}/temp_box.gro -p ${dir}/topol.top
mdrun -s ${dir}/topol.tpr -x ${dir}/minitraj.xtc -c ${dir}/res.pdb 
grep -A 4 "Energies (kJ/mol)" ${dir}/md.log > ${dir}/gromacs_ene.dat

# Convert again to phaistos

#python pdb_charmm2pdb.py ${charmm_pdb} ${phaistos_pdb}

#/Users/sandro/Programs/Phaistos/phaistos/modules/force_field/test/test_ff ${file}
#~/software/phaistos_1/phaistos/build/bin/evaluate_observable --pdb-file ${ph} \
#    --observable-opls-bond-stretch --observable-opls-angle-bend --observable-opls-torsion \
#    --observable-opls-vdw --observable-opls-charge >  ${dir}/phaistos_ene.dat

#python gromacs2human.py gromacs_ene.dat phaistos_ene.dat

