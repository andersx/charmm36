#!/bin/bash

# Remove old stuff 
rm \#*

pdb=$1
name=$2
# Create gro and top
pdb2gmx -f ${pdb} -o ${name}.gro -p ${name}.top -ignh -ter <<EOF
8
6
0
0
EOF

# preprocess
#editconf -bt cubic -f ${name}.gro -o ${name}_box.gro -c -d 1.0  
grompp -f mini_1.mdp -c ${name}.gro -p ${name}.top  -o ${name}_1.tpr
mdrun -s ${name}_1.tpr -x minitraj_1.xtc -c ${name}_gromacs.pdb 
# Creat new gro
rm \#*
pdb2gmx -f ${name}_gromacs.pdb -o ${name}_gromacs.gro -p ${name}.top <<EOF
8
6
EOF

# Create index file
make_ndx -f ${name}_gromacs.gro -o ${name}_index.ndx <<EOF
2
!2
q
EOF

grompp -f mini.mdp -c ${name}_gromacs.gro -p ${name}.top -n ${name}_index.ndx -o ${name}.tpr
mdrun -s ${name}.tpr -x minitraj_2.xtc -c ${name}_gromacs_post.pdb -e ${name}_energy.edr
g_energy -f ${name}_energy.edr
#grep -A 4 "Energies (kJ/mol)" ${dir}/md.log > ${dir}/gromacs_ene.dat
