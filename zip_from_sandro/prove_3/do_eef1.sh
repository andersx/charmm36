#!/bin/bash

pdb=$1

#rm samples/*.pdb
# Run with phaistos and produce 100 samples#
#~/Programs/Phaistos/phaistos/build/bin/phaistos --config-file config_file.txt --pdb-file ${pdb} 
rm results/compare.txt
topology=~/Programs/charmm/c36a4_rev31/toppar/top_all27_prot_na.rtf
parameter=~/Programs/charmm/c36a4_rev31/toppar/par_all27_prot_na.prm
dir=~/Programs/Phaistos/phaistos/ff_scripts/prove_3/samples/*.pdb
charmm=~/Programs/charmm/c36a4_rev31/exec/osx/charmm 
for file in $dir
do
    echo ${file}
    rm \#*
    rm temp/*
    rm results/\#*
    # convert to  charmm
    python2.6 pdb_charmm.py ${file} temp/charmm_pre.pdb
    ~/Programs/wordom_0.22-rc2/bin/wordom -mono -imol temp/charmm_pre.pdb -otrj temp/traj.dcd
    $charmm topology:${topology} parameter:${parameter} pdb:temp/charmm_pre.pdb  trajectory:temp/traj.dcd < calc_energy.inp > results/energy_charmm.txt

    python2.6 pdb_charmm2pdb.py temp/charmm_post.pdb temp/phaistos_pre.pdb
    # evaluate the energy with phaistos
    ~/Programs/Phaistos/phaistos/build_debug/test/test_ff temp/phaistos_pre.pdb temp/phaistos_post.pdb > results/energy_phaistos.txt

    python2.6 picker_2.py results/energy_phaistos.txt results/energy_gromacs.xvg results/energy_charmm.txt >> results/compare.txt
    exit
done

# compare results
#python2.6 correlate_2.py results/compare.txt