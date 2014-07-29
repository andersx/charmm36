// term_bond_stretch.h --- CHARMM36/EEF1-SB bond stretch energy term
// Copyright (C) 2014 Anders S. Christensen
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

#ifndef TERM_CHARMM_BONDSTRETCH_H
#define TERM_CHARMM_BONDSTRETCH_H

#include <string>
#include "energy/energy_term.h"
#include "parsers/topology_parser.h"
#include "parameters/bond_stretch_itp.h"

namespace phaistos {


//! bondstretch energy term
class TermCharmmBondStretch: public EnergyTermCommon<TermCharmmBondStretch, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmmBondStretch, ChainFB> EnergyTermCommon;

public:

     // Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! List of all bonded pair interactions that need to be computed
     std::vector<topology::BondedPairInteraction> bonded_pair_interactions;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmmBondStretch(ChainFB *chain,
                            const Settings &settings=Settings(),
                            RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm-bond-stretch", settings, random_number_engine) {

          // Read parameters from parameter file
          std::vector<topology::BondedPairParameter> bonded_pair_parameters 
              = topology::read_bonded_pair_parameters(charmm_constants::bond_stretch_itp);

          // Generate bond stretch terms
          this->bonded_pair_interactions = topology::generate_bonded_pair_interactions(this->chain, bonded_pair_parameters);

     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmmBondStretch(const TermCharmmBondStretch &other,
                            RandomNumberEngine *random_number_engine,
                            int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {

          // Read parameters from parameter file
          std::vector<topology::BondedPairParameter> bonded_pair_parameters 
              = topology::read_bonded_pair_parameters(charmm_constants::bond_stretch_itp);

          // Generate bond stretch terms
          this->bonded_pair_interactions = topology::generate_bonded_pair_interactions(this->chain, bonded_pair_parameters);

     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return bond stretch potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double e_bond = 0.0;

          for (unsigned int i = 0; i < this->bonded_pair_interactions.size(); i++){

               topology::BondedPairInteraction pair = this->bonded_pair_interactions[i];

               const double r = ((pair.atom1)->position - (pair.atom2)->position).norm() * charmm_constants::ANGS_TO_NM;
               const double kb = pair.kb;
               const double r0 = pair.r0;

               const double dr = r - r0;
               const double e_bond_temp = 0.5 * kb * dr * dr;

               e_bond += e_bond_temp;


               if (this->settings.debug > 1) {

                   std::cout << "# CHARMM bond-stretch-term:" 

                             << " a1: " << pair.atom1
                             << " a2: " << pair.atom2

                             << " r: : " << r * charmm_constants::NM_TO_ANGS
                             << " e_stretch : " << e_bond_temp

                             << std::endl;
                             
               }
          
          }
          if (this->settings.debug > 0) {
               printf("     bond-stretch E = %15.6f kJ/mol\n", e_bond);
               printf("     bond-stretch E = %15.6f kcal/mol\n", e_bond * charmm_constants::KJ_TO_KCAL);
          }

          return e_bond * charmm_constants::KJ_TO_KCAL;
     }

};

} // End namespace phaistos

#endif
