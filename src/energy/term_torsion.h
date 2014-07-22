// term_torsion.h --- torsion angle energy term
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

#ifndef TERM_CHARMM36_TORSION_H
#define TERM_CHARMM36_TORSION_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "parsers/topology_parser.h"
#include "math.h"

namespace phaistos {

//! Torsion energy term
class TermCharmm36Torsion: public EnergyTermCommon<TermCharmm36Torsion, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36Torsion, ChainFB> EnergyTermCommon;

     std::vector<topology::TorsionInteraction> torsion_interactions;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36Torsion(ChainFB *chain,
                        const Settings &settings=Settings(),
                        RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-torsion", settings, random_number_engine) {

          std::string filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/torsion.itp";
          std::vector<topology::TorsionParameter> torsion_parameters = topology::read_torsion_parameters(filename);
          this->torsion_interactions= topology::generate_torsion_interactions(this->chain, torsion_parameters);

     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36Torsion(const TermCharmm36Torsion &other,
                        RandomNumberEngine *random_number_engine,
                        int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {

          std::string filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/torsion.itp";
          std::vector<topology::TorsionParameter> torsion_parameters = topology::read_torsion_parameters(filename);
          this->torsion_interactions= topology::generate_torsion_interactions(this->chain, torsion_parameters);

     }
     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {


        double e_torsion = 0.0;

        for (unsigned int i = 0; i < this->torsion_interactions.size(); i++) {

            topology::TorsionInteraction torsion = this->torsion_interactions[i];

            double angle = calc_dihedral((torsion.atom1)->position,
                                         (torsion.atom2)->position,
                                         (torsion.atom3)->position,
                                         (torsion.atom4)->position);

            const double e_torsion_temp = torsion.cp * std::cos(torsion.mult * angle - torsion.phi0 * charmm36_constants::DEG_TO_RAD) + torsion.cp;

            e_torsion += e_torsion_temp;

        }

        if (settings.debug > 0) {
            printf("          torsion E = %15.6f kJ/mol\n", e_torsion);
            printf("          torsion E = %15.6f kcal/mol\n", e_torsion * charmm36_constants::KJ_TO_KCAL);
        }

        return e_torsion * charmm36_constants::KJ_TO_KCAL;

     }

};

}

#endif
