// term_imptor.h --- OPLS energy: improper torsion or out-of-plane bending energy term
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

#ifndef TERM_CHARMM36_IMPTOR_H
#define TERM_CHARMM36_IMPTOR_H

#include "energy/energy_term.h"
#include "parsers/topology_parser.h"

namespace phaistos {

//! Gromacs improper torsion Waals interaction term
class TermCharmm36Imptor: public EnergyTermCommon<TermCharmm36Imptor, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36Imptor, ChainFB> EnergyTermCommon;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     std::vector<topology::ImptorInteraction> imptor_interactions;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36Imptor(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-imptor", settings, random_number_engine) {

          std::string filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/imptor.itp";
          std::vector<topology::ImptorParameter> imptor_parameters = topology::read_imptor_parameters(filename);
          this->imptor_interactions = topology::generate_imptor_interactions(this->chain, imptor_parameters);

     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36Imptor(const TermCharmm36Imptor &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {

          std::string filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/imptor.itp";
          std::vector<topology::ImptorParameter> imptor_parameters = topology::read_imptor_parameters(filename);
          this->imptor_interactions = topology::generate_imptor_interactions(this->chain, imptor_parameters);
     }

     //! Evaluate
     //! \param move_info object containing information about last move
     //! \return improper torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energy_imptor = 0.0;

          for (unsigned int i = 0; i < this->imptor_interactions.size(); i++) {

               topology::ImptorInteraction imptor = this->imptor_interactions[i];

               const double phi = calc_dihedral((imptor.atom1)->position,
                                                (imptor.atom2)->position,
                                                (imptor.atom3)->position,
                                                (imptor.atom4)->position);

               const double dphi = phi - imptor.phi0 * charmm36_constants::DEG_TO_RAD;
               const double energy_imptor_temp = 0.5 * imptor.cp * dphi * dphi;

               energy_imptor += energy_imptor_temp;
          }

          printf("           imptor E = %15.6f kJ/mol\n", energy_imptor);
          printf("           imptor E = %15.6f kcal/mol\n", energy_imptor * charmm36_constants::KJ_TO_KCAL);

          return energy_imptor * charmm36_constants::KJ_TO_KCAL;
     }

};

}

#endif
