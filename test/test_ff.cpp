// test_opls.cpp --- Test of OPLS energy class
// Copyright (C) 2008 Kristoffer Enøe Johansson, Wouter Boomsma
//
// This file is part of Phaistos 
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.

#include "protein/chain_fb.h"
#include "energy/energy.h"

#include "energy/term_bond_stretch.h"
#include "energy/term_angle_bend.h"
#include "energy/term_torsion.h"
#include "energy/term_vdw.h"
#include "energy/term_coulomb.h"
#include "energy/term_eef1.h"


#include "energy/observable.h"
#include "energy/observable_collection.h"

using namespace std;
using namespace phaistos;
using namespace phaistos::definitions;

//! Method to evaluate a PDB file using energy terms
void test_terms(ChainFB *chain, int n=1) {


     std::cout << "Setting up force field ... " << std::endl;
     // Create Energy class 
     Energy<ChainFB> energy(chain);

     // Add terms

     energy.add_term( new TermGromacsBondStretch(chain) );
     energy.add_term( new TermGromacsAngleBend(chain) );
     energy.add_term( new TermGromacsTorsion(chain) );
     energy.add_term( new TermGromacsVdw(chain) );
     // energy.add_term( new TermOplsImptor(chain) );
     // energy.add_term( new TermGromacsCoulomb(chain) );
     // energy.add_term( new TermEef1(chain) );



     std::cout << "Init ok! " << std::endl;
     // Evaluate energy
     // for (unsigned int i = 0; i < 100; i++)
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
