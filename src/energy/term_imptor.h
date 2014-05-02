// opls_imptor.h --- OPLS energy: improper torsion or out-of-plane bending energy term
// Copyright (C) 2009-2011 Kristoffer Enøe Johansson, Wouter Boomsma
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

#ifndef TERM_GROMACS_IMPTOR_H
#define TERM_GROMACS_IMPTOR_H

#include "energy/energy_term.h"

namespace phaistos {


//! OPLS out-of-plane bending for trivalent atoms (improper torsion) energy term
class TermOplsImptor: public EnergyTermCommon<TermOplsImptor, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermOplsImptor, ChainFB> EnergyTermCommon;

     //! Number of interactions calculated
     int counter;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermOplsImptor(ChainFB *chain, const Settings &settings=Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "opls-improper-torsion", settings, random_number_engine) {

     }

     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermOplsImptor(const TermOplsImptor &other,
                    RandomNumberEngine *random_number_engine,
                    int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            counter(other.counter) {

     }

     //! Evaluate imptor energy for a atom
     //! \param atom3 central atom
     //! \return improper torsional or out-of-plane bending energy
     inline double calc_imptor_energy(Atom *atom3) {
          //     std::cout << atom3->covalent_neighbours.size() << " " << parameters.get_bond_n(atom3) << "\n";
          //     assert(atom3->covalent_neighbours.size() == 3);

          // Import protein definitions (such as residue names)
          using namespace definitions;

          // OPLS uses potonated histidines in which case there is a imptor at ND1 and
          // NE2. Phaistos has depotonated hisidines and no imptor at ND1
          // Should probably check all titratable hydrogens before continuing in case of native initialization
          // with hydrogens
          // if(atom3->atom_type == ND1 && atom3->residue->residue_type == HIS)
          if(atom3->atom_type == NE2 && atom3->residue->residue_type == HIS)
               return 0.0;

          CovalentBondIterator<ChainFB> it1(atom3, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
          int i=0;
          Atom *neighbour[3];
          for (; !it1.end(); ++it1) {
               neighbour[i] = &*it1;
               i += 1;
          }

          std::pair<ImptorParameters::Parameter, std::vector<Atom*> > ogly =
               parameters.get_param_and_sorted_atoms(neighbour[0],neighbour[1],atom3,neighbour[2]);
          ImptorParameters::Parameter param = ogly.first;
          std::vector<Atom*> sorted_atoms = ogly.second;
          double angle = calc_dihedral(sorted_atoms[0]->position,sorted_atoms[1]->position,
                                       atom3->position,sorted_atoms[3]->position);
          double e = calc_spring_energy(angle,param);
          int p1 = parameters.get_param_id(sorted_atoms[0]);
          int p2 = parameters.get_param_id(sorted_atoms[1]);
          /* int p3 = parameters.get_param_id(atom3); */
          int p4 = parameters.get_param_id(sorted_atoms[3]);

          if (p1 == p2 && p2 == p4) {
               // six-fold symmetri
               double energies[6];
               energies[0] = e;
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,e/6.0); */
               angle = calc_dihedral(sorted_atoms[1]->position,sorted_atoms[0]->position,
                                     sorted_atoms[2]->position,sorted_atoms[3]->position);
               energies[1] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,energies[1]/6.0); */
               angle = calc_dihedral(sorted_atoms[1]->position,sorted_atoms[3]->position,
                                     sorted_atoms[2]->position,sorted_atoms[0]->position);
               energies[2] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,energies[2]/6.0); */
               angle = calc_dihedral(sorted_atoms[3]->position,sorted_atoms[1]->position,
                                     sorted_atoms[2]->position,sorted_atoms[0]->position);
               energies[3] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,energies[3]/6.0); */
               angle = calc_dihedral(sorted_atoms[3]->position,sorted_atoms[0]->position,
                                     sorted_atoms[2]->position,sorted_atoms[1]->position);
               energies[4] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,energies[4]/6.0); */
               angle = calc_dihedral(sorted_atoms[0]->position,sorted_atoms[3]->position,
                                     sorted_atoms[2]->position,sorted_atoms[1]->position);
               energies[5] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,energies[5]/6.0); */
               for (i=1; i<6; i++)
                    e += energies[i];
               e /= 6.0;
          } else if (p1 == p2) {
               double energies[2];
               energies[0] = e;
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,e/2.0); */
               angle = calc_dihedral(sorted_atoms[1]->position,sorted_atoms[0]->position,
                                     sorted_atoms[2]->position,sorted_atoms[3]->position);
               energies[1] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p2,p1,p3,p4,angle*180.0/3.14159265,energies[1]/2.0); */
               e = (energies[0]+energies[1])/2.0;
          } else if (p2 == p4) {
               double energies[2];
               energies[0] = e;
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,e/2.0); */
               angle = calc_dihedral(sorted_atoms[0]->position,sorted_atoms[3]->position,
                                     sorted_atoms[2]->position,sorted_atoms[1]->position);
               energies[1] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p4,p3,p2,angle*180.0/3.14159265,energies[1]/2.0); */
               e = (energies[0]+energies[1])/2.0;
          } else if (p1 == p4) {
               double energies[2];
               energies[0] = e;
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,e/2.0); */
               angle = calc_dihedral(sorted_atoms[3]->position,sorted_atoms[1]->position,
                                     sorted_atoms[2]->position,sorted_atoms[0]->position);
               energies[1] = calc_spring_energy(angle,param);
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p4,p2,p3,p1,angle*180.0/3.14159265,energies[1]/2.0); */
               e = (energies[0]+energies[1])/2.0;
          } else {
               /* printf("%4d%4d%4d%4d %8.3f %8.5f\n",p1,p2,p3,p4,angle*180.0/3.14159265,e); */
          }
          return e;
     };

     //! Evaluate a single imptor term
     //! \param angle Dihedral angle
     //! \param param Parameter container
     inline double calc_spring_energy(double angle, ImptorParameters::Parameter param) {
          double energy=0.0;
          counter++;
          for (unsigned int i=0; i<param.size(); i++) {
               if (param.o[i] <= 0)
                    break;
               //cos is symmetric around angle so we need not flip sign when param is mirrored
               energy += param.amp[i]*( 1 + cos(param.o[i]*angle - param.ang[i]) );
          }
          return energy;
     }

     //! Evaluate
     //! \param move_info object containing information about last move
     //! \return improper torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {
          counter = 0;
          double energy_sum = 0.0;
          // all atoms
          int size = (this->chain)->size();
          for (int r = 0; r < size; r++) {
               Residue *res = &(*(this->chain))[r];
               int res_size = res->size();
               for (int a = 0; a < res_size; a++) {
                    Atom *atom = res->atoms[a];
                    // Atom must have excatly 3 neighbours (sp2 hybridized) to have
                    //  a improper torsion
                    if(parameters.get_bond_n(atom)==3)
                         energy_sum += calc_imptor_energy(atom);
               }
          }
          /* std::cout<<"Improper torsion: "<<energy_sum<<" kcal/mol  "<<counter<<" interactions\n"; */
          return energy_sum;
     };
};

}

#endif
