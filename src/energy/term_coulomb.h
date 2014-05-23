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


namespace phaistos {


//! partial coulomb interaction term
class TermCharmm36Coulomb: public EnergyTermCommon<TermCharmm36Coulomb, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36Coulomb, ChainFB> EnergyTermCommon;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     std::vector<topology::NonBondedPair> non_bonded_pairs;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36Coulomb(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-coulomb", settings, random_number_engine) {

              std::string non_bonded_filename
                  = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/vdw.itp";

              std::vector<topology::NonBondedParameter> non_bonded_parameters
                  = topology::read_nonbonded_parameters(non_bonded_filename);

              std::string non_bonded_14_filename
                  = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/vdw14.itp";

              std::vector<topology::NonBonded14Parameter> non_bonded_14_parameters =
                  topology::read_nonbonded_14_parameters(non_bonded_14_filename);

              non_bonded_pairs = topology::generate_non_bonded_pairs_cached(this->chain,
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
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            non_bonded_pairs(other.non_bonded_pairs) {}



     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double coul_energy = 0.0;
          double coul14_energy = 0.0;

          for (unsigned int i = 0; i < non_bonded_pairs.size(); i++) {

              topology::NonBondedPair pair = non_bonded_pairs[i];

               const double r_sq = ((pair.atom1)->position - (pair.atom2)->position).norm_squared();

               // This is the energy if no distance dependent di-electric constant is used.
               // const double inv_r_sq = 100.0 / (r_sq); // 100.0 due to shift to nanometers
               // const double coul_energy_temp = pair.qq * sqrt(inv_r_sq);

               // Here a distance dependent dieelectric constant of eps_r = 1.5 * r is used.
               // The factor of 10.0 here is because the pair.qq assumes distances in nanometers,
               // while the factor of 1.5 in eps_r assumes that r is in angstrom.
               const double coul_energy_temp = pair.qq / (r_sq * 1.5) * 10.0;

               if (pair.is_14_interaction) {
                    coul14_energy += coul_energy_temp;
// printf("ASC: L14COL: %5d %5d %5d  r = %8.4f  XYZ = %8.4f %8.4f %8.4f   XYZ2 = %8.4f %8.4f %8.4f   q1 = %6.3f  q2 = %6.3f  qq = %7.3f  ecoul = %7.3f\n",
// 0,0,0,
// ((pair.atom1)->position - (pair.atom2)->position).norm(),
// (pair.atom1)->position[0],
// (pair.atom1)->position[1],
// (pair.atom1)->position[2],
// (pair.atom2)->position[0],
// (pair.atom2)->position[1],
// (pair.atom2)->position[2],
// pair.q1, pair.q2, pair.qq, coul_energy_temp);

               } else {
                    coul_energy += coul_energy_temp;
// printf("ASC: LnCOL: %5d %5d %5d  r = %8.4f  XYZ = %8.4f %8.4f %8.4f   XYZ2 = %8.4f %8.4f %8.4f   q1 = %6.3f  q2 = %6.3f  qq = %7.3f  ecoul = %7.3f\n",
// 0,0,0,
// ((pair.atom1)->position - (pair.atom2)->position).norm(),
// (pair.atom1)->position[0],
// (pair.atom1)->position[1],
// (pair.atom1)->position[2],
// (pair.atom2)->position[0],
// (pair.atom2)->position[1],
// (pair.atom2)->position[2],
// pair.q1, pair.q2, pair.qq, coul_energy_temp);
               }
          }

          const double total_energy = (coul14_energy + coul_energy) / 4.184;

          printf("          Coul-14 E = %15.6f kJ/mol\n", coul14_energy);
          printf("          Coul-14 E = %15.6f kcal/mol\n", coul14_energy/4.184);
          printf("          Coul-SR E = %15.6f kJ/mol\n", coul_energy);
          printf("          Coul-SR E = %15.6f kcal/mol\n", coul_energy/4.184);
          printf("       Coul-total E = %15.6f kJ/mol\n", total_energy * 4.184);
          printf("       Coul-total E = %15.6f kcal/mol\n", total_energy);

          return total_energy;

     }
};


}

#endif
