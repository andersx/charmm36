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

#include <string>

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "charmm22_parser.h"

namespace phaistos {

//! Gromacs van der Waals interaction term
class TermGromacsVdw: public EnergyTermCommon<TermGromacsVdw, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermGromacsVdw, ChainFB> EnergyTermCommon;

     //! Number of interactions calculated
     int counter;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     std::vector<NonBondedParameter> non_bonded_parameters;
     std::vector<NonBonded14Parameter> non_bonded_14_parameters;
     std::vector<NonBondedPair> non_bonded_pairs;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermGromacsVdw(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "gromacs-vdw", settings, random_number_engine) {

              std::string non_bonded_filename = "/home/andersx/phaistos_dev/modules/gromacs/src/energy/charmm22_cmap/charmm22_vdw.itp";
              non_bonded_parameters = read_nonbonded_parameters(non_bonded_filename);

              std::string non_bonded_14_filename = "/home/andersx/phaistos_dev/modules/gromacs/src/energy/charmm22_cmap/charmm22_vdw14.itp";
              non_bonded_14_parameters = read_nonbonded_14_parameters(non_bonded_14_filename);

              // std::cout << non_bonded_parameters[1].atom_type << std::endl;
              // std::cout << non_bonded_14_parameters[1].atom_type1 << std::endl;

              non_bonded_pairs = generate_non_bonded_pairs(this->chain,
                                                           non_bonded_parameters,
                                                           non_bonded_14_parameters);

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

          double vdw_energy = 0.0;
          double coul_energy = 0.0;

          // std::cout << non_bonded_pairs.size() << " terms total" << std::endl;
          for (unsigned int i = 0; i < non_bonded_pairs.size(); i++) {

               NonBondedPair pair = non_bonded_pairs[i];

               // const double dx = (pair.atom1)->position[0] - (pair.atom2)->position[0];
               // const double dy = (pair.atom1)->position[1] - (pair.atom2)->position[1];
               // const double dz = (pair.atom1)->position[2] - (pair.atom2)->position[2];
               // const double r_sq = dx*dx + dy*dy + dz*dz;

               const double r_sq = ((pair.atom1)->position - (pair.atom2)->position).norm_squared();
               const double inv_r_sq = 100.0 / (r_sq); //shift to nanometers
               const double inv_r_sq6 = inv_r_sq * inv_r_sq * inv_r_sq;
               const double inv_r_sq12 = inv_r_sq6 * inv_r_sq6;
               const double vdw_energy_temp = pair.c12 * inv_r_sq12 - pair.c6 * inv_r_sq6;

               const double coul_energy_temp = pair.qq * sqrt(inv_r_sq);

               coul_energy += coul_energy_temp;
               vdw_energy += vdw_energy_temp;

               printf("ASC: LxCOL: %5d %5d %5d  r = %8.4f  XYZ = %8.4f %8.4f %8.4f   XYZ2 = %8.4f %8.4f %8.4f   q1 = %6.3f  q2 = %6.3f  qq = %7.3f  ecoul = %7.3f\n",
                    pair.i1, pair.i2, pair.i2 , sqrt(r_sq),
                     (pair.atom1)->position[0],  (pair.atom1)->position[1],  (pair.atom1)->position[2],
                     (pair.atom2)->position[0],  (pair.atom2)->position[1],  (pair.atom2)->position[2],
                    pair.q1, pair.q2, pair.qq, coul_energy_temp);
                  //std::cout << "ASC: " << pair.i1
                  // << "   " << pair.i2
                  // // << "   " << pair.atom1
                  // // << "   " << pair.atom2
                  // << "  q1 = " << pair.q1
                  // << "  q2 = " << pair.q2
                  // << "  qq = " << pair.qq
                  // // << "  s1 = " << pair.sigma1
                  // // << "  s2 = " << pair.sigma2
                  // // << "  e1 = " << pair.epsilon1
                  // // << "  e2 = " << pair.epsilon2
                  // << "   r = " << 10.0 / sqrt(inv_r_sq)
                  // // << std::setprecision(9) << "   c6 = "   << pair.c6
                  // // << "   c12 = "  << pair.c12
                  // // << "   evdw = "  << vdw_energy_temp
                  // << "   ecoul = "  << coul_energy_temp
                  // // << "   TOTAL = "  << vdw_energy
                  // << std::endl;

          }


          // for (AtomIterator<ChainFB, definitions::ALL> it1(*this->chain);
          //      !it1.end(); ++it1) {

          //      Atom *atom1 = &*it1;
          //      // Residue *res1 = &*it1;

          //      GromacsAtomParameter gromacs_atom_parameter1(atom1);
          //      // std::cout << atom1 << std::endl;


          // }




          // this->counter = 0;

          // std::cout << "start";
          // // Iterate all the atom pairs on the chain
          // for (AtomIterator<ChainFB, definitions::ALL> it1(*this->chain); !it1.end(); ++it1) {
          //      for(AtomIterator<ChainFB, definitions::ALL> it2(it1+1); !it2.end(); ++it2){

          //           Atom *atom1 = &*it1;
          //           Atom *atom2 = &*it2;
          //           // energy_sum += calculate_contribution(atom1, atom2);

          //           const double dx = atom1->position[0] - atom2->position[0];
          //           const double dy = atom1->position[1] - atom2->position[1];
          //           const double dz = atom1->position[2] - atom2->position[2];
          //           const double r_sq = dx*dx + dy*dy + dz*dz;
          //           const double eps = 0.0; //ASC: Epsilon
          //           const double rmin = 1.0; //ASC: Minum energy distance
          //           const double ratio = (rmin*rmin)/r_sq;
          //           const double pow6 = ratio*ratio*ratio;
          //           return (eps*(pow6*pow6-pow6));
          //      }

          // }

          // return energy_sum;

          std::cout << " COUL_TOTAL(SR) E = " << coul_energy << std::endl;
          std::cout << "  VDW_TOTAL(SR) E = " << vdw_energy << std::endl;

          const double total_energy_in_kcal_mol = (coul_energy + vdw_energy) / 4.184;
          return total_energy_in_kcal_mol;
     }
};

}

#endif
