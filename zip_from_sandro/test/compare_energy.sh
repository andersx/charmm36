#!/bin/bash

file=$1
dir=/home/sandro/projects/charmm_ph/test
topology=~/software/c36a4_rev31/toppar/top_all27_prot_na.rtf
#/home/sandro/projects/eef1/folded/input/top_all22_prot.inp
parameter=~/software/c36a4_rev31/toppar/par_all27_prot_na.prm
#/home/sandro/projects/eef1/folded/input/par_all22_prot.inp
rm \#*
rm ${dir}/md.log
rm gromacs_ene.dat 
rm charmm_ene.dat
rm ${dir}/charmmpdb.pdb
# Generate CHARMM pdb file
charmm_pdb=${dir}/charmmpdb.pdb
gro_pdb=${dir}/gropdb.pdb
phaistos_pdb=${dir}/phapdb.pdb
charmm topology:${topology} parameter:${parameter} pdbout:${charmm_pdb} \
    pdb:\"${file}\" < ${dir}/write_pdb.inp 


# first, transform to a single-frame trajectory
catdcd -otype dcd -o ${dir}/traj.dcd  -stype pdb -s ${charmm_pdb} -pdb ${charmm_pdb}
# Calculate energy with CHARMM
charmm topology:${topology} parameter:${parameter} pdb:${charmm_pdb} \
    trajectory:${dir}/traj.dcd < ${dir}/calc_energy.inp > ${dir}/charmm_out.dat

grep -a -A 9 "ENER ENR:" ${dir}/charmm_out.dat> ${dir}/charmm_ene.dat
# Calculate ENE with GROMACS

python ${dir}/pdb_charmm2gro.py ${charmm_pdb} ${gro_pdb}

pdb2gmx -f ${gro_pdb} -o ${dir}/temp.gro <<EOF 
8
6
EOF
editconf -bt cubic -f ${dir}/temp.gro -o ${dir}/temp_box.gro -c -d 1.0  
grompp -f ${dir}/mini.mdp -c ${dir}/temp_box.gro -p ${dir}/topol.top
mdrun -s ${dir}/topol.tpr -x ${dir}/minitraj.xtc -c ${dir}/res.pdb 
grep -A 4 "Energies (kJ/mol)" ${dir}/md.log > ${dir}/gromacs_ene.dat

# Convert again to phaistos

python pdb_charmm2pdb.py ${charmm_pdb} ${phaistos_pdb}
~/software/phaistos_1/phaistos/build/bin/evaluate_observable --pdb-file ${phaistos_pdb} \
    --observable-opls-bond-stretch --observable-opls-angle-bend --observable-opls-torsion \
    --observable-opls-vdw --observable-opls-charge >  ${dir}/phaistos_ene.dat

python gromacs2human.py gromacs_ene.dat charmm_ene.dat phaistos_ene.dat

