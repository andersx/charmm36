// term_vdw.h --- Van der Waals interaction energy term
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

#ifndef TERM_CHARMM36_VDW_H
#define TERM_CHARMM36_VDW_H

#include <string>

#include <boost/tokenizer.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "parsers/topology_parser.h"

namespace phaistos {


//! Gromacs van der Waals interaction term
class TermCharmm36Vdw: public EnergyTermCommon<TermCharmm36Vdw, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36Vdw, ChainFB> EnergyTermCommon;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     std::vector<topology::NonBondedParameter> non_bonded_parameters;
     std::vector<topology::NonBonded14Parameter> non_bonded_14_parameters;
     std::vector<topology::NonBondedPair> non_bonded_pairs;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36Vdw(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-vdw", settings, random_number_engine) {

              std::string non_bonded_filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/charmm22_cmap/charmm22_vdw.itp";
              non_bonded_parameters = topology::read_nonbonded_parameters(non_bonded_filename);

              std::string non_bonded_14_filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/charmm22_cmap/charmm22_vdw14.itp";
              non_bonded_14_parameters = topology::read_nonbonded_14_parameters(non_bonded_14_filename);

              non_bonded_pairs = generate_non_bonded_pairs_cached(this->chain,
                                                           non_bonded_parameters,
                                                           non_bonded_14_parameters);
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36Vdw(const TermCharmm36Vdw &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            non_bonded_pairs(other.non_bonded_pairs) {}



     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double vdw_energy = 0.0;
          double vdw14_energy = 0.0;

          // #pragma omp parallel for reduction(+:vdw_energy,vdw14_energy) schedule(static)
          for (unsigned int i = 0; i < non_bonded_pairs.size(); i++) {

              topology::NonBondedPair pair = non_bonded_pairs[i];

               const double r_sq = ((pair.atom1)->position - (pair.atom2)->position).norm_squared();
               const double inv_r_sq = 100.0 / (r_sq); //shift to nanometers
               const double inv_r_sq6 = inv_r_sq * inv_r_sq * inv_r_sq;
               const double inv_r_sq12 = inv_r_sq6 * inv_r_sq6;
               const double vdw_energy_temp = pair.c12 * inv_r_sq12 - pair.c6 * inv_r_sq6;

               if (pair.is_14_interaction) {
                    vdw14_energy += vdw_energy_temp;
               } else {
                    vdw_energy += vdw_energy_temp;
               }
          }

          // printf("           vdW-14 E = %12.4f kJ/mol\n", vdw14_energy);
          // printf("           vdW-SR E = %12.4f kJ/mol\n", vdw_energy);

          return (vdw14_energy + vdw_energy) / 4.184;
     }
};

}

#endif
