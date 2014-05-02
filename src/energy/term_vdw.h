// term_vdw.h --- Van der Waals interaction energy term
// Copyright (C) 2009-2014 Kristoffer Enøe Johansson, Wouter Boomsma,
// Anders S. Christensen
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

#ifndef TERM_GROMACS_VDW_H
#define TERM_GROMACS_VDW_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

namespace phaistos {


//! OPLS van der Waals interaction term
class TermGromacsVdw: public EnergyTermCommon<TermGromacsVdw, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermGromacsVdw, ChainFB> EnergyTermCommon;

     //! Number of interactions calculated
     int counter;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermGromacsVdw(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "gromacs-vdw", settings, random_number_engine) {

     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermGromacsVdw(const TermGromacsVdw &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            counter(other.counter) {
     }

     //! Evaluate van der Waal interaction between 2 atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \return van der Waals energy in kcal/mol
     double calculate_contribution(Atom *atom1, Atom *atom2) {

          int d = chain_distance<ChainFB>(atom1,atom2); //6s @ 500 iterations
          const double e14fac = 1.0;

          if (d > 3) {
               return calc_vdw_energy(atom1,atom2);
          } else if (d < 3) {
               return 0.0;
          } else {//d==3
               return e14fac*calc_vdw_energy(atom1,atom2);
          }
     }

     //! Evaluate Lennard-Jones potential
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param param van der Waals parameters
     //! \return lennard Jones contribution in kcal/mol
     double calc_vdw_energy(Atom *atom1, Atom *atom2) {
          counter++;
          const double dx = atom1->position[0] - atom2->position[0];
          const double dy = atom1->position[1] - atom2->position[1];
          const double dz = atom1->position[2] - atom2->position[2];
          const double r_sq = dx*dx + dy*dy + dz*dz;
          const double eps = 0.0; //ASC: Epsilon
          const double rmin = 1.0; //ASC: Minum energy distance
          const double ratio = (rmin*rmin)/r_sq;
          const double pow6 = ratio*ratio*ratio;
          return (eps*(pow6*pow6-pow6));
     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energy_sum = 0.0;
          this->counter = 0;

          // Iterate all the atom pairs on the chain
          for (AtomIterator<ChainFB, definitions::ALL> it1(*this->chain); !it1.end(); ++it1) {
               for(AtomIterator<ChainFB, definitions::ALL> it2(it1+1); !it2.end(); ++it2){

                    Atom *atom1 = &*it1;
                    Atom *atom2 = &*it2;
                    energy_sum += calculate_contribution(atom1, atom2);
               }

          }

          return energy_sum;
     }
};

}

#endif
