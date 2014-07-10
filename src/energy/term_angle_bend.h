// term_angle_bend.h --- angle-bend energy term
// Copyright (C) 2014 Anders S. Christensn
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

#ifndef TERM_CHARMM36_ANGLEBEND_H
#define TERM_CHARMM36_ANGLEBEND_H

#include <string>

#include <boost/tokenizer.hpp>

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "parsers/topology_parser.h"

namespace phaistos {


//! OPLS anglebend energy term - base class containing all functionality
class TermCharmm36AngleBend: public EnergyTermCommon<TermCharmm36AngleBend, ChainFB> {

protected:

    //! For convenience, define local EnergyTermCommon
    typedef phaistos::EnergyTermCommon<TermCharmm36AngleBend, ChainFB> EnergyTermCommon;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     std::vector<topology::AngleBendInteraction> angle_bend_interactions;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36AngleBend(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-angle-bend", settings, random_number_engine) {

          std::string filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/angle_bend.itp";
          std::vector<topology::AngleBendParameter> angle_bend_parameters = topology::read_angle_bend_parameters(filename);
          this->angle_bend_interactions = topology::generate_angle_bend_interactions(this->chain, angle_bend_parameters);
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36AngleBend(const TermCharmm36AngleBend &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {

          std::string filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/angle_bend.itp";
          std::vector<topology::AngleBendParameter> angle_bend_parameters = topology::read_angle_bend_parameters(filename);
          this->angle_bend_interactions = topology::generate_angle_bend_interactions(this->chain, angle_bend_parameters);

     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return angle bend potential energy of the chain in the object
     double evaluate(MoveInfo *move_info = NULL) {

          double energy_sum = 0.0;
          double energy_angle = 0.0;
          double energy_urey = 0.0;

          for (unsigned int i = 0; i < this->angle_bend_interactions.size(); i++){

               topology::AngleBendInteraction interaction = this->angle_bend_interactions[i];

               const double theta = calc_angle((interaction.atom1)->position,
                                               (interaction.atom2)->position,
                                               (interaction.atom3)->position);

               // Angle bend part
               const double dtheta = theta - interaction.theta0 * M_PI / 180.0;
               const double energy_angle_bend_temp = 0.5 * interaction.k0 * dtheta * dtheta; 

               energy_angle += energy_angle_bend_temp;

               // Urey-Bradley part
               const double r13 = ((interaction.atom1)->position - (interaction.atom3)->position).norm() * charmm36_constants::ANGS_TO_NM;
               const double dr = r13 - interaction.r13;
               const double energy_urey_bradley_temp = 0.5 * interaction.kub * dr * dr;

               energy_urey += energy_urey_bradley_temp;
               energy_sum += energy_angle_bend_temp + energy_urey_bradley_temp;

          }

          if (settings.debug >= 0) {

               printf("       angle-bend E = %15.6f kJ/mol\n", energy_angle);
               printf("       angle-bend E = %15.6f kcal/mol\n", 
                       energy_angle * charmm36_constants::KJ_TO_KCAL);
               printf("     urey-bradley E = %15.6f kJ/mol\n", energy_urey);
               printf("     urey-bradley E = %15.6f kcal/mol\n",
                       energy_urey * charmm36_constants::KJ_TO_KCAL);
               printf(" angle-bend-total E = %15.6f kJ/mol\n", energy_sum);
               printf(" angle-bend-total E = %15.6f kcal/mol\n",
                       energy_sum * charmm36_constants::KJ_TO_KCAL);

          }

          return energy_sum * charmm36_constants::KJ_TO_KCAL;

     }


};

} // End namespace phaistos
#endif
