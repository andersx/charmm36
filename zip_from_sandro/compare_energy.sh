#!/bin/bash

file=$1
dir=/home/sandro/projects/charmm_ph/test
rm \#*
# Convert file
pdb2gmx -f ${file} -o ${dir}/temp.gro <<EOF
9
6
EOF
exit

editconf -bt cubic -f ${dir}/temp.gro -o ${dir}/temp_box.gro -c -d 1.0
grompp -f ${dir}/mini.mdp -c ${dir}/temp_box.gro -p ${dir}/topol.top
mdrun -s ${dir}/topol.tpr -o ${dir}/topol.trj -c ${dir}/res.pdb

