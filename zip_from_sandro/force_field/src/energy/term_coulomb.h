// term_coulomb.h ---  coulomb-coulomb interaction energy term
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

#ifndef TERM_COULOMB_H
#define TERM_COULOMB_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"


namespace phaistos {



//! partial coulomb interaction term
class TermCoulomb: public EnergyTermCommon<TermCoulomb, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCoulomb, ChainFB> EnergyTermCommon;


     //! Number of interactions in the last evaluation
     int counter;

     // float fast_inv_sqrt(float x);

public:

     //! Local settings class
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:

          //! Scaling for 1-4 interactions
          double E14FAC;

          //! Dielectric constant
          double EPS;

          //! Whether to use distance-dependent dielectric constant
          bool RDIE;

          //! Constructor
          Settings(double E14FAC=1.0,double EPS=1.0, bool RDIE=false)
               : E14FAC(E14FAC),EPS(EPS),RDIE(RDIE){}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "E14FAC-coulomb:" << settings.E14FAC << "\n";
               o << "Dielectric constant:" << settings.EPS << "\n";
               o << "Distance-dependent-dielectric:" << settings.RDIE << "\n";
               o << static_cast<const EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }
     } settings;    //!< Local settings object

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param chain Molecule chain
     TermCoulomb(ChainFB *chain,
                 const Settings &settings = Settings(),
                 RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "energy-coulomb", settings, random_number_engine) {

     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCoulomb(const TermCoulomb &other,
                    RandomNumberEngine *random_number_engine,
                    int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            counter(other.counter){

     }

     //! Evaluate coulomb interaction between 2 atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param chg1 Coulomb of atom1
     //! \param chg2 Coulomb of atom2
     //! \return Coulomb energy for atom pair
     double calculate_contribution(Atom *atom1, Atom *atom2) {

          int d = chain_distance<ChainFB>(atom1,atom2); //6s @ 500 iterations

          if (d > 3) {
               return calc_coulomb_energy(atom1,atom2);
          } else if (d < 3) {
               return 0.0;
          } else {//d==3
               return (settings.E14FAC)*calc_coulomb_energy(atom1,atom2);
          }
     }

     //! Evaluate coulomb interaction between 2 atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param chg1 Coulomb of atom1
     //! \param chg2 Coulomb of atom2
     //! \return Coulomb energy for atom pair
     double calculate_rdie_contribution(Atom *atom1, Atom *atom2) {

          int d = chain_distance<ChainFB>(atom1,atom2); //6s @ 500 iterations

          if (d > 3) {
               return calc_rdie_coulomb_energy(atom1,atom2);
          } else if (d < 3) {
               return 0.0;
          } else {//d==3
               return (settings.E14FAC)*calc_rdie_coulomb_energy(atom1,atom2);
          }
     }

     //! Evaluate a coulomb interaction between two atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param chg1 Coulomb of atom1
     //! \param chg2 Coulomb of atom2
     //! \return Coulomb energy for atom pair
     double calc_coulomb_energy(Atom *atom1, Atom *atom2) {

          counter++;
          const double chg1 = 1.0;
          const double chg2 = 1.0;
          double r = (atom1->position - atom2->position).norm();
          return (chg1*chg2/(settings.EPS*r));
     }
     //! Evaluate a coulomb interaction between two atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param chg1 Coulomb of atom1
     //! \param chg2 Coulomb of atom2
     //! \return Coulomb energy for atom pair
     double calc_rdie_coulomb_energy(Atom *atom1, Atom *atom2) {

          counter++;
          const double chg1 = 1.0;
          const double chg2 = 1.0;
          double r = (atom1->position - atom2->position).norm();
          return (chg1*chg2/(settings.EPS*r*r));
     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energy_sum=0.0;
          counter=0;

          // Duplicated code - but saves many if statements evaluations
          if(settings.RDIE==false){

               // Iterate all the atom pairs on the chain
               for (AtomIterator<ChainFB, definitions::ALL> it1(*this->chain); !it1.end(); ++it1) {
                    for(AtomIterator<ChainFB, definitions::ALL> it2(it1+1); !it2.end(); ++it2){

                         Atom *atom1 = &*it1;
                         Atom *atom2 = &*it2;
                         energy_sum += calculate_contribution(atom1, atom2);
                    }
               }

          } else {
               // Iterate all the atom pairs on the chain
               for (AtomIterator<ChainFB, definitions::ALL> it1(*this->chain); !it1.end(); ++it1) {
                    for(AtomIterator<ChainFB, definitions::ALL> it2(it1+1); !it2.end(); ++it2){

                         Atom *atom1 = &*it1;
                         Atom *atom2 = &*it2;
                         energy_sum += calculate_rdie_contribution(atom1, atom2);
                    }
               }
          }

          return energy_sum;
     }

};


}

#endif
