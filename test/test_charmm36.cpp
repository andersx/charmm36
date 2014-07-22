// test_charmm36.cpp --- Test of CHARMM36/EEF1-SB energy class
// Copyright (C) 2014 Anders Steen Christensen
//
// This file is part of PHAISTOS 
//
// PHAISTOS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// PHAISTOS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with PHAISTOS.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <string.h>

#include <boost/tokenizer.hpp>

#include "protein/chain_fb.h"
#include "protein/definitions.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "energy/energy.h"

#include "energy/observable.h"
#include "energy/observable_collection.h"

#include "energy/term_bond_stretch.h"
#include "energy/term_bonded_cached.h"
#include "energy/term_non_bonded_cached.h"
#include "energy/term_angle_bend.h"
#include "energy/term_torsion.h"
#include "energy/term_imptor.h"
#include "energy/term_vdw.h"
#include "energy/term_coulomb.h"
#include "energy/term_eef1.h"
#include "energy/term_cmap.h"
using namespace std;
using namespace phaistos;
using namespace phaistos::definitions;

//! Method to evaluate a PDB file using energy terms
void test_terms(ChainFB *chain, int n=1) {


     std::cout << "Setting up force field ... " << std::endl;
     // Create Energy class 
     Energy<ChainFB> energy(chain);

     // Add terms
     TermCharmm36BondStretch::Settings settings_bond_stretch;
     TermCharmm36AngleBend::Settings settings_angle_bend;
     TermCharmm36Torsion::Settings settings_torsion;
     TermCharmm36Vdw::Settings settings_vdw;
     TermCharmm36Imptor::Settings settings_imptor;
     TermCharmm36Coulomb::Settings settings_coulomb;
     TermCharmm36Eef1::Settings settings_eef1;
     TermCharmm36Cmap::Settings settings_cmap;

     settings_bond_stretch.debug = 1;
     settings_angle_bend.debug = 1;
     settings_torsion.debug = 1;
     settings_vdw.debug = 1;
     settings_imptor.debug = 1;
     settings_coulomb.debug = 1;
     settings_eef1.debug = 1;
     settings_cmap.debug = 1;

     energy.add_term( new TermCharmm36BondStretch(chain, settings_bond_stretch) );
     energy.add_term( new TermCharmm36AngleBend(chain, settings_angle_bend) );
     energy.add_term( new TermCharmm36Torsion(chain, settings_torsion) );
     energy.add_term( new TermCharmm36Vdw(chain, settings_vdw) );
     energy.add_term( new TermCharmm36Imptor(chain, settings_imptor) );
     energy.add_term( new TermCharmm36Coulomb(chain, settings_coulomb) );
     energy.add_term( new TermCharmm36Eef1(chain, settings_eef1) );
     energy.add_term( new TermCharmm36Cmap(chain, settings_cmap) );
     energy.add_term( new TermCharmm36BondedCached(chain) );
     energy.add_term( new TermCharmm36NonBondedCached(chain) );



     std::cout << "Init ok! " << std::endl;

     // Evaluate energy

     energy.evaluate();
     // Output
     cout << energy << endl;;

     
}



int main(int argc, char *argv[]) {

     if (argc < 2) {
          cout << "USAGE: ./test_ff <pdb-file>" <<endl;;
          exit(1);
     }

     // Create chain from PDB filename
     string pdb_filename = argv[1];
     ChainFB chain(pdb_filename, ALL_ATOMS);
     // Add atoms missing in the pdb structure
     //chain.add_atoms(ALL_PHYSICAL_ATOMS);

     test_terms(&chain, 1);
}
