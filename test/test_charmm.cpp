// test_charmm.cpp --- Test of CHARMM36/EEF1-SB energy class
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
#include "energy/term_non_bonded.h"
#include "energy/term_non_bonded_cached.h"
#include "energy/term_angle_bend.h"
#include "energy/term_torsion.h"
#include "energy/term_improper_torsion.h"
#include "energy/term_vdw.h"
#include "energy/term_coulomb.h"
#include "energy/term_implicit_solvent.h"
#include "energy/term_cmap.h"

//! Method to evaluate a PDB file using energy terms
void test_terms(phaistos::ChainFB *chain, int debug_level) {

     using namespace phaistos;
     using namespace definitions;

     std::cout << "Setting up force field ... " << std::endl;

     // Create Energy class 
     Energy<ChainFB> energy(chain);

     // Create settings object
     TermCharmmBondStretch::Settings settings_bond_stretch;
     TermCharmmAngleBend::Settings settings_angle_bend;
     TermCharmmTorsion::Settings settings_torsion;
     TermCharmmVdw::Settings settings_vdw;
     TermCharmmImproperTorsion::Settings settings_improper_torsion;
     TermCharmmCoulomb::Settings settings_coulomb;
     TermCharmmImplicitSolvent::Settings settings_implicit_solvent;
     TermCharmmCmap::Settings settings_cmap;

     // Set debug-level (which controls output)
     settings_bond_stretch.debug = debug_level;
     settings_angle_bend.debug = debug_level;
     settings_torsion.debug = debug_level;
     settings_vdw.debug = debug_level;
     settings_improper_torsion.debug = debug_level;
     settings_coulomb.debug = debug_level;
     settings_implicit_solvent.debug = debug_level;
     settings_cmap.debug = debug_level;

     // Add terms to energy object
     energy.add_term(new TermCharmmBondStretch(chain, settings_bond_stretch));
     energy.add_term(new TermCharmmAngleBend(chain, settings_angle_bend));
     energy.add_term(new TermCharmmTorsion(chain, settings_torsion));
     energy.add_term(new TermCharmmVdw(chain, settings_vdw));
     energy.add_term(new TermCharmmImproperTorsion(chain, settings_improper_torsion));
     energy.add_term(new TermCharmmCoulomb(chain, settings_coulomb));
     energy.add_term(new TermCharmmImplicitSolvent(chain, settings_implicit_solvent));
     energy.add_term(new TermCharmmCmap(chain, settings_cmap));
     energy.add_term(new TermCharmmNonBonded(chain));
     energy.add_term(new TermCharmmBondedCached(chain));
     energy.add_term(new TermCharmmNonBondedCached(chain));

     // Evaluate energy
     energy.evaluate();

     // Output
     std::cout << energy << std::endl;;

}


int main(int argc, char *argv[]) {

     using namespace phaistos;
     using namespace definitions;

     int debug_level = 1;

     if (argc == 3) {

         debug_level = atoi(argv[2]);
     }

     if (argc < 2) {
         std::cout << "USAGE: ./test_charmm <pdb-file> [debug-level (=1 by default)]" << std::endl;;
          exit(1);
     }



     // Create chain from PDB filename
     std::string pdb_filename = argv[1];
     ChainFB chain(pdb_filename, ALL_ATOMS);

     test_terms(&chain, debug_level);
}


