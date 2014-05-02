#!/bin/bash

pdb=$1

#rm samples/*.pdb
# Run with phaistos and produce 100 samples
#~/Programs/Phaistos/phaistos/build/bin/phaistos --config-file config_file.txt --pdb-file ${pdb} 
rm results/compare.txt
dir=~/Programs/Phaistos/phaistos/ff_scripts/prove_3/samples/*.pdb
for file in $dir
do
    rm \#*
    rm temp/*
    rm results/\#*
    
    # evaluate the energy with phaistos
    ~/Programs/Phaistos/phaistos/build_debug/test/test_ff ${file} temp/temp_phaistos.pdb > results/energy_phaistos.txt
    python2.6 pdb_gromacs.py temp/temp_phaistos.pdb
    pdb2gmx -f temp/temp_phaistos_mod.pdb  -o temp/temp.gro -p temp/temp.top <<EOF
8
6
EOF
    grompp -f mini.mdp -c temp/temp.gro -p temp/temp.top  -o temp/temp.tpr
    mdrun -s temp/temp.tpr -x temp/minitraj.xtc -c temp/temp_gromacs.pdb -e temp/temp_ener.edr
    g_energy -f temp/temp_ener.edr -o results/energy_gromacs.xvg <<EOF
7
9
EOF

    python2.6 picker.py results/energy_phaistos.txt results/energy_gromacs.xvg >> results/compare.txt
done

# compare results
python2.6 correlate.py results/compare.txt