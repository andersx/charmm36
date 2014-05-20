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

     std::vector<topology::AngleBendPair> angle_bend_pairs;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36AngleBend(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-angle-bend", settings, random_number_engine) {

          std::string filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/charmm22_cmap/charmm22_angle_bend.itp";
          std::vector<topology::AngleBendParameter> angle_bend_parameters = topology::read_angle_bend_parameters(filename);
          this->angle_bend_pairs = topology::generate_angle_bend_pairs(this->chain, angle_bend_parameters);
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36AngleBend(const TermCharmm36AngleBend &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            angle_bend_pairs(other.angle_bend_pairs) {}

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return angle bend potential energy of the chain in the object
     double evaluate(MoveInfo *move_info = NULL) {

          double energy_sum = 0.0;

          for (unsigned int i = 0; i < this->angle_bend_pairs.size(); i++){

               topology::AngleBendPair pair = this->angle_bend_pairs[i];

               const double theta = calc_angle((pair.atom1)->position,
                                               (pair.atom2)->position,
                                               (pair.atom3)->position);

               // Angle bend part
               const double dtheta = theta - pair.theta0 * M_PI / 180.0;
               const double energy_angle_bend_temp = 0.5 * pair.k0 * dtheta * dtheta; 

               // Urey-Bradley part
               const double r13 = ((pair.atom1)->position - (pair.atom3)->position).norm() / 10.0;
               const double dr = r13 - pair.r13;
               const double energy_urey_bradley_temp = 0.5 * pair.kub * dr * dr;

               energy_sum += energy_angle_bend_temp + energy_urey_bradley_temp;

          }

          // printf("     urey-bradley E = %12.4f kJ/mol\n", energy_sum);

          return energy_sum / 4.184;

     }


};

} // End namespace phaistos
#endif
