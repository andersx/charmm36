// term_coulomb.h ---  coulomb-coulomb interaction energy term
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

#ifndef TERM_CHARMM36_COULOMB_H
#define TERM_CHARMM36_COULOMB_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"

#include "parameters/vdw14_itp.h"
#include "parameters/vdw_itp.h"

namespace phaistos {


//! partial coulomb interaction term
class TermCharmm36Coulomb: public EnergyTermCommon<TermCharmm36Coulomb, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36Coulomb, ChainFB> EnergyTermCommon;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! List that holds all the interactions and parameters that need to be computed:w
     std::vector<topology::NonBondedInteraction> non_bonded_interactions;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36Coulomb(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-coulomb", settings, random_number_engine) {

              std::vector<topology::NonBondedParameter> non_bonded_parameters
                  = topology::read_nonbonded_parameters(charmm36_constants::vdw_itp);

              std::vector<topology::NonBonded14Parameter> non_bonded_14_parameters =
                  topology::read_nonbonded_14_parameters(charmm36_constants::vdw14_itp);

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
     TermCharmm36Coulomb(const TermCharmm36Coulomb &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {


              std::vector<topology::NonBondedParameter> non_bonded_parameters
                  = topology::read_nonbonded_parameters(charmm36_constants::vdw_itp);

              std::vector<topology::NonBonded14Parameter> non_bonded_14_parameters =
                  topology::read_nonbonded_14_parameters(charmm36_constants::vdw14_itp);

              this->non_bonded_interactions = 
                  topology::generate_non_bonded_interactions(this->chain,
                                                             non_bonded_parameters,
                                                             non_bonded_14_parameters);
     }



     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double coul_energy = 0.0;
          double coul14_energy = 0.0;

          for (unsigned int i = 0; i < this->non_bonded_interactions.size(); i++) {

              topology::NonBondedInteraction interaction = this->non_bonded_interactions[i];

               const double r_sq = ((interaction.atom1)->position - (interaction.atom2)->position).norm_squared();

               // This is the energy if no distance dependent di-electric constant is used.
               // const double inv_r_sq = 100.0 / (r_sq); // 100.0 due to shift to nanometers
               // const double coul_energy_temp = pair.qq * sqrt(inv_r_sq);

               // Here a distance dependent dieelectric constant of eps_r = 1.5 * r is used.
               // The factor of 10.0 here is because the pair.qq assumes distances in nanometers,
               // while the factor of 1.5 in eps_r assumes that r is in angstrom, so only one r in r^2 
               // must be converted.
               const double coul_energy_temp = interaction.qq / (r_sq * 1.5) * charmm36_constants::NM_TO_ANGS;

               if (interaction.is_14_interaction) {

                   coul14_energy += coul_energy_temp;

                   if (this->settings.debug > 1) {

                       std::cout << "# CHARMM36 coulomb-14:";
                  
                    }
               } else {
 
                   coul_energy += coul_energy_temp;

                    if (this->settings.debug > 1) {

                       std::cout << "# CHARMM36 coulomb:";
                  
                    }
               }

               if (this->settings.debug > 1) {

                       std::cout << " a1: " << interaction.atom1
                                 << " a2: " << interaction.atom2

                                 << " q1*q2: " <<  interaction.qq

                                 << " r: " << std::sqrt(r_sq)

                                 << " e_coul: " << coul_energy_temp

                                 << std::endl;
                }
          }

          const double total_energy = (coul14_energy + coul_energy) * charmm36_constants::KJ_TO_KCAL;

          if (this->settings.debug > 0) {
               printf("          Coul-14 E = %15.6f kJ/mol\n", coul14_energy);
               printf("          Coul-14 E = %15.6f kcal/mol\n", coul14_energy * charmm36_constants::KJ_TO_KCAL);
               printf("          Coul-SR E = %15.6f kJ/mol\n", coul_energy);
               printf("          Coul-SR E = %15.6f kcal/mol\n", coul_energy * charmm36_constants::KJ_TO_KCAL);
               printf("       Coul-total E = %15.6f kJ/mol\n", total_energy * charmm36_constants::KCAL_TO_KJ);
               printf("       Coul-total E = %15.6f kcal/mol\n", total_energy);
          }
          return total_energy;

     }
};


}

#endif
