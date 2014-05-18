// term_bond_stretch.h --- bond-stretch energy term
// Copyright (C) 2009-2011 Kristoffer Enøe Johansson
//
// This file is part of Phaistos
//
// Phaistos is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Phaistos is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef TERM_CHARMM36_BONDSTRETCH_H
#define TERM_CHARMM36_BONDSTRETCH_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "parsers/topology_parser.h"

namespace phaistos {


//! bondstretch energy term
class TermCharmm36BondStretch: public EnergyTermCommon<TermCharmm36BondStretch, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36BondStretch, ChainFB> EnergyTermCommon;

public:

     // Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     std::vector<topology::BondedPair> bonded_pairs;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36BondStretch(ChainFB *chain,
                            const Settings &settings=Settings(),
                            RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-bond-stretch", settings, random_number_engine) {

          std::string filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/charmm22_cmap/charmm22_bond.itp";

          std::vector<topology::BondedPairParameter> bonded_pair_parameters 
              = topology::read_bonded_pair_parameters(filename);

          this->bonded_pairs = topology::generate_bonded_pairs(this->chain, bonded_pair_parameters);

     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36BondStretch(const TermCharmm36BondStretch &other,
                            RandomNumberEngine *random_number_engine,
                            int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
          bonded_pairs(other.bonded_pairs) {}

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return bond stretch potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double e_bond = 0.0;

          for (unsigned int i = 0; i < this->bonded_pairs.size(); i++){

               topology::BondedPair pair = this->bonded_pairs[i];

               const double r = ((pair.atom1)->position - (pair.atom2)->position).norm() / 10.0;
               const double kb = pair.kb;
               const double r0 = pair.r0;

               const double dr = r - r0;
               const double e_bond_temp = 0.5 * kb * dr * dr;

               e_bond += e_bond_temp;

               // printf("ASC: BONDSTRETCH    dr = %14.10f   r0 = %14.10f   kb = %14.10f  vbond = %14.10f\n", dr, r0, kb, e_bond_temp);
          }

          // printf("     bond-stretch E = %12.4f kJ/mol\n", e_bond);

          return e_bond / 4.184;
     }

};

}

#endif
