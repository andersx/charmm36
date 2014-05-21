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
                  = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/charmm22_cmap/charmm22_vdw.itp";

              std::vector<topology::NonBondedParameter> non_bonded_parameters
                  = topology::read_nonbonded_parameters(non_bonded_filename);

              std::string non_bonded_14_filename
                  = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/charmm22_cmap/charmm22_vdw14.itp";

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
               const double inv_r_sq = 100.0 / (r_sq); //shift to nanometers
               const double coul_energy_temp = pair.qq * sqrt(inv_r_sq);

               if (pair.is_14_interaction) {
                    coul14_energy += coul_energy_temp;
               } else {
                    coul_energy += coul_energy_temp;
               }
          }

          printf("          Coul-14 E = %12.4f kJ/mol\n", coul14_energy);
          printf("          Coul-SR E = %12.4f kJ/mol\n", coul_energy);

          return (coul14_energy + coul_energy) / 4.184;

     }
};


}

#endif
