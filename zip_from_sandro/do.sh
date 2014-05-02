#!/bin/bash

pdb=$1

rm samples/*.pdb
# Run with phaistos and produce 100 samples
~/Programs/Phaistos/phaistos/build/bin/phaistos --config-file config_file.txt --pdb-file ${pdb} 

dir=~/Programs/Phaistos/phaistos/ff_scripts/prove_3/${system}/samples/*.pdb
for file in $dir
do
    rm \#*
    rm temp/\#*
    
    echo 8 6 | pdb2gmx -f ${file} -o temp/temp.gro -p temp/temp.top -ignh 
    grompp -f mini.mdp -c temp/temp.gro -p temp/temp.top  -o temp/temp.tpr
    mdrun -s temp/temp.tpr -x temp/minitraj.xtc -c $ temp/temp_gromacs.pdb 

    # evaluate the energy with phaistos
    ~/Programs/Phaistos/phaistos/build_debug/test_ff temp/temp_gromacs.pdb >> results/energy_${system}.txt

	
    echo 8 6 | pdb2gmx -f temp/temp_phaistos.pdb -o temp/temp_phaistos.gro -p temp/temp.top 
    grompp -f mini.mdp -c temp/temp_phaistos.gro -p temp/temp.top  -o temp/temp.tpr
    mdrun -s temp/temp.tpr -x minitraj.xtc -c $ temp_gromacs.pdb 
    g_energy -f ${name}_energy.edr
    exit
done

# compare results