// term_vdw.h --- CHARMM36/EEF1-SB van der Waals interaction energy term
// Copyright (C) 2014 Anders S. Christensen
//
// This file is part of PHAISTOS
//
// PHAISTOS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// PHAISTOS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Phaistos.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef TERM_CHARMM_VDW_H
#define TERM_CHARMM_VDW_H

#include <string>

#include <boost/tokenizer.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"

#include "parsers/topology_parser.h"

#include "parameters/vdw14_itp.h"
#include "parameters/vdw_itp.h"

namespace phaistos {


//! Gromacs van der Waals interaction term
class TermCharmmVdw: public EnergyTermCommon<TermCharmmVdw, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmmVdw, ChainFB> EnergyTermCommon;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     std::vector<topology::NonBondedInteraction> non_bonded_interactions;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmmVdw(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm-vdw", settings, random_number_engine) {

              std::vector<topology::NonBondedParameter> non_bonded_parameters
                  = topology::read_nonbonded_parameters(charmm_constants::vdw_itp);

              std::vector<topology::NonBonded14Parameter> non_bonded_14_parameters =
                  topology::read_nonbonded_14_parameters(charmm_constants::vdw14_itp);

              this->non_bonded_interactions = 
                  topology::generate_non_bonded_interactions(this->chain,
                                                             non_bonded_parameters,
                                                             non_bonded_14_parameters);
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmmVdw(const TermCharmmVdw &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {


              std::vector<topology::NonBondedParameter> non_bonded_parameters
                  = topology::read_nonbonded_parameters(charmm_constants::vdw_itp);

              std::vector<topology::NonBonded14Parameter> non_bonded_14_parameters =
                  topology::read_nonbonded_14_parameters(charmm_constants::vdw14_itp);

              this->non_bonded_interactions = 
                  topology::generate_non_bonded_interactions(this->chain,
                                                             non_bonded_parameters,
                                                             non_bonded_14_parameters);
     }



     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double vdw_energy = 0.0;
          double vdw14_energy = 0.0;

          for (unsigned int i = 0; i < this->non_bonded_interactions.size(); i++) {

              topology::NonBondedInteraction interaction = this->non_bonded_interactions[i];

               const double r2 = ((interaction.atom1)->position - (interaction.atom2)->position).norm_squared();
               const double inv_r2 = charmm_constants::NM2_TO_ANGS2 / r2;
               const double inv_r6 = inv_r2 * inv_r2 * inv_r2;
               const double inv_r12 = inv_r6 * inv_r6;
               const double vdw_energy_temp = interaction.c12 * inv_r12 - interaction.c6 * inv_r6;

               if (interaction.is_14_interaction) {

                    vdw14_energy += vdw_energy_temp;

                    if (this->settings.debug > 1) {

                       std::cout << "# CHARMM vdw-14:";
                  
                    }
 
               } else {
 
                   vdw_energy += vdw_energy_temp;

                    if (this->settings.debug > 1) {

                       std::cout << "# CHARMM vdw:";
                  
                    }
               }

               if (this->settings.debug > 1) {

                       std::cout << " a1: " << interaction.atom1
                                 << " a2: " << interaction.atom2

                                 << " c6: " <<  interaction.c6
                                 << " c12: " <<  interaction.c12

                                 << " r: " << std::sqrt(r2)

                                 << " e_vdw: " << vdw_energy_temp

                                 << std::endl;
                }



          }

          const double total_energy = (vdw14_energy + vdw_energy) * charmm_constants::KJ_TO_KCAL;

          if (this->settings.debug > 0) {
              printf("           vdW-14 E = %15.6f kJ/mol\n", vdw14_energy);
              printf("           vdW-14 E = %15.6f kcal/mol\n", vdw14_energy * charmm_constants::KJ_TO_KCAL);
              printf("           vdW-SR E = %15.6f kJ/mol\n", vdw_energy);
              printf("           vdW-SR E = %15.6f kcal/mol\n", vdw_energy * charmm_constants::KJ_TO_KCAL);
              printf("        vdW-total E = %15.6f kJ/mol\n", total_energy * charmm_constants::KCAL_TO_KJ);
              printf("        vdW-total E = %15.6f kcal/mol\n", total_energy);
          }

          return total_energy;

     }
};

}

#endif
