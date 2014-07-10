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

     //! Struct that holds lists of all interaction the needs 
     //! to be computed if a residue is changed.
     struct BondedCachedResidue {

          std::vector<topology::AngleBendInteraction> angle_bend_interactions;
          std::vector<topology::BondedPairInteraction> bonded_pair_interactions;
          std::vector<topology::ImptorInteraction> imptor_interactions;
          std::vector<topology::TorsionInteraction> torsion_interactions;
          topology::CmapInteraction cmap_interaction;

          // If the CMAP correction has to be computed for this residue, 
          // i.e. if it is not first or last residue in the chain.
          bool has_cmap;

          // Current energy of the chain.
          double energy_new;

          // Backup of last energy, which is used in caching.
          double energy_old;

     };

     //! Table which contains the CMAP correction tables
     //! in Gromacs' format
     std::vector<std::vector<double> > cmap_data;

     //! Vector containing a list of all interaction 
     std::vector<BondedCachedResidue> bonded_cached_residues;

     double energy_new;
     double energy_old;

     int start_index;
     int end_index;


     //! Setup 
     void setup_caches() {

          // Get CMAP data from the Gromacs code.
          this->cmap_data = charmm36_cmap::setup_cmap();
          std::vector<topology::CmapInteraction> cmap_interactions 
              = topology::generate_cmap_interactions(this->chain);

          std::string filename;

          // Read angle-bend parameters.
          filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/angle_bend.itp";
          std::vector<topology::AngleBendParameter> angle_bend_parameters 
              = topology::read_angle_bend_parameters(filename);

          std::vector<topology::AngleBendInteraction> angle_bend_interactions 
              = topology::generate_angle_bend_interactions(this->chain, angle_bend_parameters);

          // Read bond-stretch parameters.
          filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/bond_stretch.itp";
          std::vector<topology::BondedPairParameter> bonded_pair_parameters
              = topology::read_bonded_pair_parameters(filename);

          std::vector<topology::BondedPairInteraction> bonded_pair_interactions
              = topology::generate_bonded_pair_interactions(this->chain, bonded_pair_parameters);

          // Read improper torsion parameters.
          filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/imptor.itp";
          std::vector<topology::ImptorParameter> imptor_parameters 
              = topology::read_imptor_parameters(filename);
          std::vector<topology::ImptorInteraction> imptor_interactions
              = topology::generate_imptor_interactions(this->chain, imptor_parameters);

          // Get proper torsion parameters.
          filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/torsion.itp";
          std::vector<topology::TorsionParameter> torsion_parameters
              = topology::read_torsion_parameters(filename);

          std::vector<topology::TorsionInteraction> torsion_interactions
              = topology::generate_torsion_interactions(this->chain, torsion_parameters);

          // Make room in cache vector for all residues.
          this->bonded_cached_residues.resize(this->chain->size());

          // Sort angle_bend_pairs
          for (unsigned int i = 0; i < angle_bend_interactions.size(); i++){

               topology::AngleBendInteraction interaction = angle_bend_interactions[i];

               int index = std::min((interaction.atom1)->residue->index,
                                    (interaction.atom2)->residue->index);

               index = std::min((interaction.atom3)->residue->index,
                                index);

               this->bonded_cached_residues[index].angle_bend_interactions.push_back(interaction);

          }

          // Sort bonded_pairs
          for (unsigned int i = 0; i < bonded_pair_interactions.size(); i++){

               topology::BondedPairInteraction interaction = bonded_pair_interactions[i];

               int index = std::min((interaction.atom1)->residue->index,
                                    (interaction.atom2)->residue->index);

               this->bonded_cached_residues[index].bonded_pair_interactions.push_back(interaction);

          }

          // Sort improper torsions
          for (unsigned int i = 0; i < imptor_interactions.size(); i++){

               topology::ImptorInteraction interaction = imptor_interactions[i];

               int index = std::min((interaction.atom1)->residue->index,
                                    (interaction.atom2)->residue->index);

               index = std::min((interaction.atom3)->residue->index,
                                index);

               index = std::min((interaction.atom4)->residue->index,
                                index);

               this->bonded_cached_residues[index].imptor_interactions.push_back(interaction);

          }

          // Sort proper torsions
          for (unsigned int i = 0; i < torsion_interactions.size(); i++){

               topology::TorsionInteraction interaction = torsion_interactions[i];

               int index = std::min((interaction.atom1)->residue->index,
                                    (interaction.atom2)->residue->index);

               index = std::min((interaction.atom3)->residue->index,
                                index);

               index = std::min((interaction.atom4)->residue->index,
                                index);

               this->bonded_cached_residues[index].torsion_interactions.push_back(interaction);

          }

          for (unsigned int i = 0; i < this->bonded_cached_residues.size(); i ++) {
               this->bonded_cached_residues[i].has_cmap = false;
          }

          // Sort CMAP parameters
          for (unsigned int i = 0; i < cmap_interactions.size(); i++){

               topology::CmapInteraction interaction = cmap_interactions[i];

               int index = interaction.residue_index;
               this->bonded_cached_residues[index].cmap_interaction = interaction;
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

     }



     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36BondedCached(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-bonded-cached", settings, random_number_engine) {


         setup_caches();
     }

     double calculate_cached_residue_energy(BondedCachedResidue &cached_residue) {


          double energy_sum = 0.0;

          // Sort angle_bend_pairs
          for (unsigned int i = 0; i < cached_residue.angle_bend_interactions.size(); i++){

               topology::AngleBendInteraction interaction = cached_residue.angle_bend_interactions[i];

               const double theta = calc_angle((interaction.atom1)->position,
                                               (interaction.atom2)->position,
                                               (interaction.atom3)->position);

               // Angle bend part
               const double dtheta = theta - interaction.theta0 * charmm36_constants::DEG_TO_RAD;
               const double energy_angle_bend_temp = 0.5 * interaction.k0 * dtheta * dtheta;

               // Urey-Bradley part
               const double r13 = ((interaction.atom1)->position - (interaction.atom3)->position).norm() * charmm36_constants::ANGS_TO_NM;
               const double dr = r13 - interaction.r13;
               const double energy_urey_bradley_temp = 0.5 * interaction.kub * dr * dr;

               energy_sum += energy_angle_bend_temp + energy_urey_bradley_temp;

          }

          for (unsigned int i = 0; i < cached_residue.bonded_pair_interactions.size(); i++){

               topology::BondedPairInteraction interaction = cached_residue.bonded_pair_interactions[i];

               const double r = ((interaction.atom1)->position - (interaction.atom2)->position).norm() * charmm36_constants::ANGS_TO_NM;
               const double kb = interaction.kb;
               const double r0 = interaction.r0;

               const double dr = r - r0;
               const double e_bond_temp = 0.5 * kb * dr * dr;

               energy_sum += e_bond_temp;


          }

          for (unsigned int i = 0; i < cached_residue.imptor_interactions.size(); i++){

               topology::ImptorInteraction interaction = cached_residue.imptor_interactions[i];

               const double phi = calc_dihedral((interaction.atom1)->position,
                                                (interaction.atom2)->position,
                                                (interaction.atom3)->position,
                                                (interaction.atom4)->position);

               const double dphi = phi - interaction.phi0 * charmm36_constants::DEG_TO_RAD;
               const double energy_imptor_temp = 0.5 * interaction.cp * dphi * dphi;

               energy_sum += energy_imptor_temp;


          }

          for (unsigned int i = 0; i < cached_residue.torsion_interactions.size(); i++){

               topology::TorsionInteraction interaction = cached_residue.torsion_interactions[i];

               double angle = calc_dihedral((interaction.atom1)->position,
                                            (interaction.atom2)->position,
                                            (interaction.atom3)->position,
                                            (interaction.atom4)->position);

               const double e_torsion_temp = interaction.cp * std::cos(interaction.mult * angle - interaction.phi0 * charmm36_constants::DEG_TO_RAD) + interaction.cp;

               energy_sum += e_torsion_temp;

          }

          if (cached_residue.has_cmap) {
                const int residue_index = cached_residue.cmap_interaction.residue_index;
                const unsigned int cmap_type_index = cached_residue.cmap_interaction.cmap_type_index;

                const double phi = (*(this->chain))[residue_index].get_phi();
                const double psi = (*(this->chain))[residue_index].get_psi();

                energy_sum += charmm36_cmap::cmap_energy(phi, psi, cmap_type_index, this->cmap_data);
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
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {

          setup_caches();
     }

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

          return this->energy_new * charmm36_constants::KJ_TO_KCAL;
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
