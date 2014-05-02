// term_bond_stretch.h --- bond-stretch energy term
// Copyright (C) 2009-2011 Kristoffer Enøe Johansson
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

#ifndef TERM_BONDSTRETCH_H
#define TERM_BONDSTRETCH_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"

namespace phaistos {


//! bondstretch energy term
class TermBondStretch: public EnergyTermCommon<TermBondStretch, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermBondStretch, ChainFB> EnergyTermCommon;

     //! Number of interactions calculated
     int counter;

public:

     // Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermBondStretch(ChainFB *chain,
                         const Settings &settings=Settings(),
                         RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "bond-stretch", settings, random_number_engine) {

     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermBondStretch(const TermBondStretch &other,
                         RandomNumberEngine *random_number_engine,
                         int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            counter(other.counter) {
     }

     //! Evaluate anglebend energy for a bond (atom1-atom2)
     //! \param atom2 Second atom
     //! \return Bond stretch energy to previous neighbours
     inline double calc_bondstretch_energy(Atom *atom2) {
          double energy = 0.0;
          double len = 1.5;
          double k = 1.0;
          CovalentBondIterator<ChainFB> it1(atom2, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
          for (; !it1.end(); ++it1) {
               Atom *atom1 = &*it1;

               // Only calculate the distance to covalently bonded that are before atom2 in the chain
               // (otherwise all distance would be considered twice)
               if (atom1->residue->index < atom2->residue->index ||
                   (atom1->residue->index == atom2->residue->index && atom1->index < atom2->index)) {
                    double length = (atom1->position - atom2->position).norm();
                    energy += calc_spring_energy(length,len,k);
               }
          }

          return energy;
     }

     //! Evaluate the potetial energy in a spring
     //! \param x Distance
     //! \param x_eq Equilibrium distance
     //! \param k Spring constant
     //! \return Spring energy
     inline double calc_spring_energy(double x, double x_eq, double k) {
          counter++;
          double dx = x - x_eq;
          return (k*dx*dx);
     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return bond stretch potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energySum=0.0;
          counter=0;

          //+1 to avoid getNeighbout(-1) error from N[0]
          AtomIterator<ChainFB,definitions::ALL> it(*(this->chain)); ++it;
          for (; !it.end(); ++it) {
               Atom &atom = *it;
               energySum += calc_bondstretch_energy(&atom);
          }
          return energySum;
     }

};

}

#endif
