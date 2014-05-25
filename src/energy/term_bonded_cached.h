// term_bonded_cached.h -- All bonded terms -- cached version.
// Copyright (C) 2014 Anders S. Christensen
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

#ifndef TERM_CHARMM36_BONDED_CACHED_H
#define TERM_CHARMM36_BONDED_CACHED_H

#include <string>
#include <limits>

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "parsers/topology_parser.h"

#include "kernels/kernel_cmap.h"
#include "term_cmap_tables.h"

namespace phaistos {

//! CHARMM36 bonded terms -- cached.
class TermCharmm36BondedCached: public EnergyTermCommon<TermCharmm36BondedCached, ChainFB> {

protected:

    //! For convenience, define local EnergyTermCommon
    typedef phaistos::EnergyTermCommon<TermCharmm36BondedCached, ChainFB> EnergyTermCommon;

public:


     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     struct BondedCachedResidue {

          std::vector<topology::AngleBendPair> angle_bend_pairs;
          std::vector<topology::BondedPair> bonded_pairs;
          std::vector<topology::Imptor> imptors;
          std::vector<topology::DihedralAngleType9> dihedral_angles;
          topology::CmapPair cmap_pair;

          bool has_cmap;

          double energy_new;
          double energy_old;

     };

     //! Table which contains the CMAP correction tables
     //! in Gromacs' 
     std::vector<std::vector<double> > cmapdata;

     //! Vector containing all terms in the CMAP correction.

     std::vector<BondedCachedResidue> bonded_cached_residues;

     double energy_new;
     double energy_old;

     int start_index;
     int end_index;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36BondedCached(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-bonded-cached", settings, random_number_engine) {

          // Get CMAP data from the Gromacs code.
          this->cmapdata = charmm36_cmap::setup_cmap();
          std::vector<topology::CmapPair> cmap_pairs = topology::generate_cmap_pairs(this->chain);

          // Read angle-bend parameters.
          std::string filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/angle_bend.itp";
          std::vector<topology::AngleBendParameter> angle_bend_parameters = topology::read_angle_bend_parameters(filename);

          std::vector<topology::AngleBendPair> angle_bend_pairs = topology::generate_angle_bend_pairs(this->chain, angle_bend_parameters);

          // Read bond-stretch parameters.
          filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/bond_stretch.itp";

          std::vector<topology::BondedPairParameter> bonded_pair_parameters
              = topology::read_bonded_pair_parameters(filename);

          std::vector<topology::BondedPair> bonded_pairs = topology::generate_bonded_pairs(this->chain, bonded_pair_parameters);

          // Read improper torsion parameters.
          filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/imptor.itp";
          std::vector<topology::DihedralType2Parameter> dihedral_type_2_parameters = topology::read_dihedral_type_2_parameters(filename);
          std::vector<topology::Imptor> imptors = topology::generate_imptors(this->chain, dihedral_type_2_parameters);

          // Get proper torsion parameters.
          filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/torsion.itp";
          std::vector<topology::DihedralType9Parameter> dihedral_type_9_parameters = topology::read_dihedral_type_9_parameters(filename);

          std::vector<topology::DihedralAngleType9> dihedral_angles = topology::generate_dihedral_pairs(this->chain, dihedral_type_9_parameters);

          // Make room in cache vector for all residues.
          this->bonded_cached_residues.resize(this->chain->size());

          // Sort angle_bend_pairs
          for (unsigned int i = 0; i < angle_bend_pairs.size(); i++){

               topology::AngleBendPair pair = angle_bend_pairs[i];

               int index = std::numeric_limits<int>::max();

               if ((pair.atom1)->residue->index < index) index = (pair.atom1)->residue->index;
               if ((pair.atom2)->residue->index < index) index = (pair.atom2)->residue->index;
               if ((pair.atom3)->residue->index < index) index = (pair.atom3)->residue->index;

               this->bonded_cached_residues[index].angle_bend_pairs.push_back(pair);

          }

          // Sort bonded_pairs
          for (unsigned int i = 0; i < bonded_pairs.size(); i++){

               topology::BondedPair pair = bonded_pairs[i];

               int index = std::numeric_limits<int>::max();

               if ((pair.atom1)->residue->index < index) index = (pair.atom1)->residue->index;
               if ((pair.atom2)->residue->index < index) index = (pair.atom2)->residue->index;

               this->bonded_cached_residues[index].bonded_pairs.push_back(pair);

          }

          // Sort improper torsions
          for (unsigned int i = 0; i < imptors.size(); i++){

               topology::Imptor pair = imptors[i];

               int index = std::numeric_limits<int>::max();

               if ((pair.atom1)->residue->index < index) index = (pair.atom1)->residue->index;
               if ((pair.atom2)->residue->index < index) index = (pair.atom2)->residue->index;
               if ((pair.atom3)->residue->index < index) index = (pair.atom3)->residue->index;
               if ((pair.atom4)->residue->index < index) index = (pair.atom4)->residue->index;

               this->bonded_cached_residues[index].imptors.push_back(pair);

          }

          // Sort proper torsions
          for (unsigned int i = 0; i < dihedral_angles.size(); i++){

               topology::DihedralAngleType9 pair = dihedral_angles[i];

               int index = std::numeric_limits<int>::max();

               if ((pair.atom1)->residue->index < index) index = (pair.atom1)->residue->index;
               if ((pair.atom2)->residue->index < index) index = (pair.atom2)->residue->index;
               if ((pair.atom3)->residue->index < index) index = (pair.atom3)->residue->index;
               if ((pair.atom4)->residue->index < index) index = (pair.atom4)->residue->index;

               this->bonded_cached_residues[index].dihedral_angles.push_back(pair);

         }

          for (unsigned int i = 0; i < this->bonded_cached_residues.size(); i ++) {
               this->bonded_cached_residues[i].has_cmap = false;
          }
          // Sort CMAP parameters
          for (unsigned int i = 0; i < cmap_pairs.size(); i++){

               topology::CmapPair pair = cmap_pairs[i];

               int index = pair.residue_index;
               this->bonded_cached_residues[index].cmap_pair = pair;
               this->bonded_cached_residues[index].has_cmap = true;

         }

         this->energy_new  = 0.0;
         this->energy_old  = 0.0;

         for (unsigned int i = 0; i < this->bonded_cached_residues.size(); i ++) {

             double residue_energy = calculate_cached_residue_energy(bonded_cached_residues[i]);
             this->energy_new  += residue_energy;
             this->energy_old  += residue_energy;

             this->bonded_cached_residues[i].energy_new = residue_energy;
             this->bonded_cached_residues[i].energy_old = residue_energy;

         }

        // std::cout << this->energy_new << std::endl;
        // std::cout << this->energy_new / 4.184 << std::endl;

     }

     double calculate_cached_residue_energy(BondedCachedResidue &cached_residue) {


          double energy_sum = 0.0;

          // Sort angle_bend_pairs
          for (unsigned int i = 0; i < cached_residue.angle_bend_pairs.size(); i++){

               topology::AngleBendPair pair = cached_residue.angle_bend_pairs[i];

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

          // Sort bonded_pairs
          for (unsigned int i = 0; i < cached_residue.bonded_pairs.size(); i++){

               topology::BondedPair pair = cached_residue.bonded_pairs[i];

               const double r = ((pair.atom1)->position - (pair.atom2)->position).norm() / 10.0;
               const double kb = pair.kb;
               const double r0 = pair.r0;

               const double dr = r - r0;
               const double e_bond_temp = 0.5 * kb * dr * dr;

               energy_sum += e_bond_temp;


          }

          // Sort improper torsions
          for (unsigned int i = 0; i < cached_residue.imptors.size(); i++){

               topology::Imptor imptor = cached_residue.imptors[i];

               const double phi = calc_dihedral((imptor.atom1)->position,
                                                (imptor.atom2)->position,
                                                (imptor.atom3)->position,
                                                (imptor.atom4)->position);

               const double dphi = phi - imptor.phi0 / 180.0 * M_PI;
               const double energy_imptor_temp = 0.5 * imptor.cp * dphi * dphi;

               energy_sum += energy_imptor_temp;


          }

          // Sort proper torsions
          for (unsigned int i = 0; i < cached_residue.dihedral_angles.size(); i++){

               topology::DihedralAngleType9 dihedral = cached_residue.dihedral_angles[i];

               double angle = calc_dihedral((dihedral.atom1)->position,
                                            (dihedral.atom2)->position,
                                            (dihedral.atom3)->position,
                                            (dihedral.atom4)->position);

               const double e_torsion_temp = dihedral.cp * cos(dihedral.mult * angle - dihedral.phi0 / 180.0 * M_PI) + dihedral.cp;

               energy_sum += e_torsion_temp;

          }

          if (cached_residue.has_cmap) {
                const int residue_index = cached_residue.cmap_pair.residue_index;
                const unsigned int cmap_type_index = cached_residue.cmap_pair.cmap_type_index;
                // const unsigned int cmap_type_index = 0;

                const double phi = (*(this->chain))[residue_index].get_phi();
                const double psi = (*(this->chain))[residue_index].get_psi();

                energy_sum += kernel::cmap_energy(phi, psi, cmap_type_index, this->cmapdata);
          }

          return energy_sum;


     }


     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36BondedCached(const TermCharmm36BondedCached &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            cmapdata(other.cmapdata),
            bonded_cached_residues(other.bonded_cached_residues),
            energy_new(other.energy_new),
            energy_old(other.energy_old),
            start_index(other.start_index),
            end_index(other.end_index) {}

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return angle bend potential energy of the chain in the object
     double evaluate(MoveInfo *move_info = NULL) {


          // Indexes of first and last residue for which the
          // position of atoms have changed since last move.
          this->start_index = 0;
          this->end_index = this->chain->size() - 1;

          // If these are set explicitly by the move, read these here.
          if (move_info) {
              this->start_index = std::max(0, move_info->modified_positions_start - 1);
              this->end_index = std::min(move_info->modified_positions_end + 1, this->chain->size() - 1);
          }

          // Local delta energy required for OpenMP -- can't just write to this->energy_new.
          double delta_energy_local = 0.0;

          #pragma omp parallel for reduction(+:delta_energy_local) schedule(static)
          for (int i = this->start_index; i < this->end_index+1; i ++) {

               const double residue_energy = 
                    calculate_cached_residue_energy(bonded_cached_residues[i]);

               this->bonded_cached_residues[i].energy_new = residue_energy;

               delta_energy_local += this->bonded_cached_residues[i].energy_new
                                   - this->bonded_cached_residues[i].energy_old;
          }

          this->energy_new  += delta_energy_local;

          // double energy_check = 0.0;
          // for (unsigned int i = 0; i < this->bonded_cached_residues.size(); i ++) {
          //      energy_check += calculate_cached_residue_energy(bonded_cached_residues[i]);
          // }
          // printf("ASC: Delta check - actual:  %24.14f\n", energy_check - this->energy_new);

          return this->energy_new/4.184;
     }


     void accept() {

          // If move is accepted, backup energies in all pairs that were recomputed
          for (int i = this->start_index; i < this->end_index+1; i ++) {

               this->bonded_cached_residues[i].energy_old
                   = this->bonded_cached_residues[i].energy_new;
          }

          this->energy_old = this->energy_new;
    }


    void reject() {

        // If move is accepted, restore energies in all pairs that were recomputed
          for (int i = this->start_index; i < this->end_index+1; i ++) {

               this->bonded_cached_residues[i].energy_new
                   = this->bonded_cached_residues[i].energy_old;
          }
          this->energy_new = this->energy_old;
    }


};

} // End namespace phaistos
#endif
