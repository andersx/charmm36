// term_improper_torsion.h --- Improper torsion or out-of-plane bending energy term
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

#ifndef TERM_CHARMM36_IMPROPER_TORSION_H
#define TERM_CHARMM36_IMPROPER_TORSION_H

#include "energy/energy_term.h"
#include "parsers/topology_parser.h"
#include "parameters/imptor_itp.h"

namespace phaistos {

//! Gromacs improper torsion Waals interaction term
class TermCharmm36ImproperTorsion: public EnergyTermCommon<TermCharmm36ImproperTorsion, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36ImproperTorsion, ChainFB> EnergyTermCommon;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! List of interactions that need to be computed
     std::vector<topology::ImproperTorsionInteraction> improper_torsion_interactions;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36ImproperTorsion(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-improper-torsion", settings, random_number_engine) {

          std::vector<topology::ImproperTorsionParameter> improper_torsion_parameters 
                    = topology::read_improper_torsion_parameters(charmm36_constants::imptor_itp);

          this->improper_torsion_interactions 
                    = topology::generate_improper_torsion_interactions(this->chain, improper_torsion_parameters);

     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36ImproperTorsion(const TermCharmm36ImproperTorsion &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {

          std::vector<topology::ImproperTorsionParameter> improper_torsion_parameters 
                    = topology::read_improper_torsion_parameters(charmm36_constants::imptor_itp);

          this->improper_torsion_interactions 
                    = topology::generate_improper_torsion_interactions(this->chain, improper_torsion_parameters);

     }

     //! Evaluate
     //! \param move_info object containing information about last move
     //! \return improper torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energy_improper_torsion = 0.0;

          for (unsigned int i = 0; i < this->improper_torsion_interactions.size(); i++) {

               topology::ImproperTorsionInteraction interaction = this->improper_torsion_interactions[i];

               const double phi = calc_dihedral((interaction.atom1)->position,
                                                (interaction.atom2)->position,
                                                (interaction.atom3)->position,
                                                (interaction.atom4)->position);

               const double dphi = phi - interaction.phi0 * charmm36_constants::DEG_TO_RAD;
               const double energy_improper_torsion_temp = 0.5 * interaction.cp * dphi * dphi;

               energy_improper_torsion += energy_improper_torsion_temp;

               if (this->settings.debug > 1) {

                   std::cout << "# CHARMM36 improper-torsion:" 

                             << " a1: " << interaction.atom1
                             << " a2: " << interaction.atom2
                             << " a3: " << interaction.atom3
                             << " a4: " << interaction.atom4

                             << " angle: " << phi * charmm36_constants::RAD_TO_DEG
                             << " e_improper_torsion: " <<  energy_improper_torsion_temp

                             << std::endl;
                }
          }

          if (this->settings.debug > 0) {

               printf(" improper-torsion E = %15.6f kJ/mol\n", energy_improper_torsion);
               printf(" improper-torsion E = %15.6f kcal/mol\n", energy_improper_torsion * charmm36_constants::KJ_TO_KCAL);
          }

          return energy_improper_torsion * charmm36_constants::KJ_TO_KCAL;
     }

};

}

#endif
