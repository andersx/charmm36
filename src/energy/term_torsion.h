// term_torsion.h --- torsion angle energy term
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

#ifndef TERM_GROMACS_TORSION_H
#define TERM_GROMACS_TORSION_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"

namespace phaistos {

//! Torsion energy term
class TermGromacsTorsion: public EnergyTermCommon<TermGromacsTorsion, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermGromacsTorsion, ChainFB> EnergyTermCommon;

     //! Number of interactions calculated
     int counter;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermGromacsTorsion(ChainFB *chain,
                        const Settings &settings=Settings(),
                        RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "gromacs-torsion", settings, random_number_engine) {

     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermGromacsTorsion(const TermGromacsTorsion &other,
                        RandomNumberEngine *random_number_engine,
                        int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            counter(other.counter) {

     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          return 0;

     }





//     //! Evaluate torsion energy for a bond (atom2-atom3)
//     //! \param atom2 First atom defining the bond
//     //! \param atom3 Second atom defining the bond
//     //! \return Torsional energy of that bond
//     inline double calc_torsion_energy(Atom *atom2, Atom *atom3) {
//          double angle,energy = 0.0;
//          CovalentBondIterator<ChainFB> it1(atom2, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
//          for (; !it1.end(); ++it1) {
//               Atom *atom1 = &*it1;
//               if (atom1 == atom3)
//                    continue;
//               CovalentBondIterator<ChainFB> it4(atom3, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
//               for (; !it4.end(); ++it4) {
//                    Atom *atom4 = &*it4;
//                    if (atom4 == atom2)
//                         continue;
//                    double param = 1.0;
//                    angle = calc_dihedral(atom1->position,atom2->position,
//                                          atom3->position,atom4->position);
//                    energy += calc_spring_energy(angle,param);
//               }
//          }
//          return energy;
//     }
//
//     //! Evaluate a single torsional term
//     //! \param angle Dihedral angle
//     //! \param param Parameter container
//     inline double calc_spring_energy(double angle, double param) {
//          double energy=0.0;
//          counter++;
//          //cos is symmetric around angle so we need not flip sign when param is mirrored
//          energy += param*( 1 + cos(param*angle - param) );
//
//          return energy;
//     }
//
//     //! Evaluate chain energy
//     //! \param move_info object containing information about last move
//     //! \return torsional potential energy of the chain in the object
//     double evaluate(MoveInfo *move_info=NULL) {
//
//          double energy_sum=0.0;
//          counter=0;
//
//          // all torsions
//          int size = (this->chain)->size();
//          for (int r=0; r<size; r++) {
//               Residue *res = &(*(this->chain))[r];
//               int res_size = res->size();
//               for (int a=0; a<res_size; a++) {
//                    Atom *atom2 = res->atoms[a];
//                    CovalentBondIterator<ChainFB> it(atom2, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
//                    for (; !it.end(); ++it) {
//                         Atom *atom3 = &*it;
//                         if(atom2->index < atom3->index)
//                              energy_sum += calc_torsion_energy(atom2,atom3);
//                    }
//               }
//          }
//
//          /* // phaistos degrees of freedom only */
//          /* DofIterator::angle_selection_enum dofs = DofIterator::DIHEDRAL_DOFS + */
//          /*      DofIterator::CHI_ANGLES + DofIterator::N_DIHEDRAL; */
//          /* DofIterator dofIt((*chain)(0,N),DofIterator::DIHEDRAL,dofs); */
//          /* DofIterator end = chain->dofIteratorEnd(); */
//          /* for (; dofIt!=end; ++dofIt) { */
//          /*      Atom *atom3 = dofIt.getAtom(); */
//          /*      Atom *atom1init,*atom2,*atom4init; */
//          /*      atom3->get_dihedral_atoms(&atom1init,&atom2,&atom3,&atom4init); */
//          /*      energy_sum += calc_torsion_energy(atom2,atom3); */
//          /* } */
//
//          //energy_sum *= parameters.torsion_unit;
//
//          /* std::cout<<"Torsional Angle "<<energy_sum<<" kcal/mol "<<counter<<" interactions\n"; */
//
//          return energy_sum;
//     }
};

}

#endif
