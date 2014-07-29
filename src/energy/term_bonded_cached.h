// term_bonded_cached.h -- CHARMM36/EEF1-SB cached version of all bonded terms.
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

#ifndef TERM_CHARMM_BONDED_CACHED_H
#define TERM_CHARMM_BONDED_CACHED_H

#include <string>

#include "energy/energy_term.h"
#include "parsers/topology_parser.h"
#include "term_cmap_tables.h"
#include "parameters/angle_bend_itp.h"
#include "parameters/bond_stretch_itp.h"
#include "parameters/imptor_itp.h"
#include "parameters/torsion_itp.h"


namespace phaistos {

//! CHARMM bonded terms -- cached.
class TermCharmmBondedCached: public EnergyTermCommon<TermCharmmBondedCached, ChainFB> {

protected:

    //! For convenience, define local EnergyTermCommon
    typedef phaistos::EnergyTermCommon<TermCharmmBondedCached, ChainFB> EnergyTermCommon;

public:

     //! Local settings class
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:

          // Flags to ignore individual terms
          bool ignore_bond_angles;
          bool ignore_bond_stretch;
          bool ignore_torsion_angles;
          bool ignore_improper_torsion_angles;
          bool ignore_cmap_correction;

          //! Constructor
          Settings(bool ignore_bond_angles=false,
                   bool ignore_bond_stretch=false,
                   bool ignore_torsion_angles=false,
                   bool ignore_improper_torsion_angles=false,
                   bool ignore_cmap_correction=false)
               : ignore_bond_angles(ignore_bond_angles),
                 ignore_bond_stretch(ignore_bond_stretch),
                 ignore_torsion_angles(ignore_torsion_angles),
                 ignore_improper_torsion_angles(ignore_improper_torsion_angles),
                 ignore_cmap_correction(ignore_cmap_correction) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "ignore-bond-angles:" << settings.ignore_bond_angles << "\n";
               o << "ignore-bond-stretch:" << settings.ignore_bond_stretch << "\n";
               o << "ignore-torsion-angles:" << settings.ignore_torsion_angles << "\n";
               o << "ignore-improper-torsion-angles:" << settings.ignore_improper_torsion_angles << "\n";
               o << "ignore-cmap-correction:" << settings.ignore_cmap_correction << "\n";
               o << static_cast<const EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }
     } settings;    //!< Local settings object

     //! Struct that holds lists of all interaction the needs 
     //! to be computed if a residue is changed.
     struct BondedCachedResidue {

          std::vector<topology::AngleBendInteraction> angle_bend_interactions;
          std::vector<topology::BondedPairInteraction> bonded_pair_interactions;
          std::vector<topology::ImproperTorsionInteraction> improper_torsion_interactions;
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

     //! Energy after the current move
     double energy_new;

     //! Backup of energy before current move
     double energy_old;

     //! Index of first residue that was moved in current move
     int start_index;

     //! Index of last residue that was moved in current move
     int end_index;

     //! Flag to keep track of none-moves
     bool none_move;

     //! Setup 
     void setup_caches() {

          // Get CMAP data from the Gromacs code.
          this->cmap_data = charmm_cmap::setup_cmap();
          std::vector<topology::CmapInteraction> cmap_interactions 
              = topology::generate_cmap_interactions(this->chain);

          // Read angle-bend parameters.
          std::vector<topology::AngleBendParameter> angle_bend_parameters = 
                    topology::read_angle_bend_parameters(charmm_constants::angle_bend_itp);

          std::vector<topology::AngleBendInteraction> angle_bend_interactions 
              = topology::generate_angle_bend_interactions(this->chain, angle_bend_parameters);

          // Read bond-stretch parameters.
          std::vector<topology::BondedPairParameter> bonded_pair_parameters 
              = topology::read_bonded_pair_parameters(charmm_constants::bond_stretch_itp);

          std::vector<topology::BondedPairInteraction> bonded_pair_interactions
              = topology::generate_bonded_pair_interactions(this->chain, bonded_pair_parameters);

          // Read improper torsion parameters.
          std::vector<topology::ImproperTorsionParameter> improper_torsion_parameters 
                    = topology::read_improper_torsion_parameters(charmm_constants::imptor_itp);

          std::vector<topology::ImproperTorsionInteraction> improper_torsion_interactions
              = topology::generate_improper_torsion_interactions(this->chain, improper_torsion_parameters);

          // Get proper torsion parameters.
          std::vector<topology::TorsionParameter> torsion_parameters 
                    = topology::read_torsion_parameters(charmm_constants::torsion_itp);

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
          for (unsigned int i = 0; i < improper_torsion_interactions.size(); i++){

               topology::ImproperTorsionInteraction interaction = improper_torsion_interactions[i];

               int index = std::min((interaction.atom1)->residue->index,
                                    (interaction.atom2)->residue->index);

               index = std::min((interaction.atom3)->residue->index,
                                index);

               index = std::min((interaction.atom4)->residue->index,
                                index);

               this->bonded_cached_residues[index].improper_torsion_interactions.push_back(interaction);

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

          // Initialize CMAP flags to false
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

         // Initialize energies
         this->energy_new  = 0.0;
         this->energy_old  = 0.0;

         // Initialize each cache
         for (unsigned int i = 0; i < this->bonded_cached_residues.size(); i ++) {

             double residue_energy = calculate_cached_residue_energy(bonded_cached_residues[i]);
             this->energy_new  += residue_energy;
             this->energy_old  += residue_energy;

             this->bonded_cached_residues[i].energy_new = residue_energy;
             this->bonded_cached_residues[i].energy_old = residue_energy;

         }

     }

     //! Calculate energy of a cached residue object
     //! \param cached_residue A residue object for which the energy is calculated
     //! \returns The energy of the residue
     double calculate_cached_residue_energy(BondedCachedResidue &cached_residue) {

          // Initialize residue energy
          double energy_sum = 0.0;

          // Calculate bond angle terms
          if (!(this->settings.ignore_bond_angles)) {

              for (unsigned int i = 0; i < cached_residue.angle_bend_interactions.size(); i++){

                   topology::AngleBendInteraction interaction = cached_residue.angle_bend_interactions[i];

                   const double theta = calc_angle((interaction.atom1)->position,
                                                   (interaction.atom2)->position,
                                                   (interaction.atom3)->position);

                   // Angle bend part
                   const double dtheta = theta - interaction.theta0 * charmm_constants::DEG_TO_RAD;
                   const double energy_angle_bend_temp = 0.5 * interaction.k0 * dtheta * dtheta;

                   // Urey-Bradley part
                   const double r13 = ((interaction.atom1)->position - (interaction.atom3)->position).norm() * charmm_constants::ANGS_TO_NM;
                   const double dr = r13 - interaction.r13;
                   const double energy_urey_bradley_temp = 0.5 * interaction.kub * dr * dr;

                   energy_sum += energy_angle_bend_temp + energy_urey_bradley_temp;
               }
          }

          // Calculate bond stretch terms
          if (!(this->settings.ignore_bond_stretch)) {

               for (unsigned int i = 0; i < cached_residue.bonded_pair_interactions.size(); i++){

                    topology::BondedPairInteraction interaction = cached_residue.bonded_pair_interactions[i];

                    const double r = ((interaction.atom1)->position - (interaction.atom2)->position).norm() * charmm_constants::ANGS_TO_NM;
                    const double kb = interaction.kb;
                    const double r0 = interaction.r0;

                    const double dr = r - r0;
                    const double e_bond_temp = 0.5 * kb * dr * dr;

                    energy_sum += e_bond_temp;

               }
          }

          // Calculate improper torsion terms
          if (!(this->settings.ignore_improper_torsion_angles)) {

               for (unsigned int i = 0; i < cached_residue.improper_torsion_interactions.size(); i++){

                    topology::ImproperTorsionInteraction interaction = cached_residue.improper_torsion_interactions[i];

                    const double phi = calc_dihedral((interaction.atom1)->position,
                                                     (interaction.atom2)->position,
                                                     (interaction.atom3)->position,
                                                     (interaction.atom4)->position);

                    const double dphi = phi - interaction.phi0 * charmm_constants::DEG_TO_RAD;
                    const double energy_improper_torsion_temp = 0.5 * interaction.cp * dphi * dphi;

                    energy_sum += energy_improper_torsion_temp;

               }
          }

          // Calculate torsion terms
          if (!(this->settings.ignore_torsion_angles)) {

               for (unsigned int i = 0; i < cached_residue.torsion_interactions.size(); i++){

                    topology::TorsionInteraction interaction = cached_residue.torsion_interactions[i];

                    double angle = calc_dihedral((interaction.atom1)->position,
                                                 (interaction.atom2)->position,
                                                 (interaction.atom3)->position,
                                                 (interaction.atom4)->position);

                    const double e_torsion_temp = interaction.cp * std::cos(interaction.mult * angle - interaction.phi0 * charmm_constants::DEG_TO_RAD) 
                                                    + interaction.cp;

                    energy_sum += e_torsion_temp;
              }

          }

          // Calculate CMAP correction terms
          if (!(this->settings.ignore_cmap_correction)) {
               if (cached_residue.has_cmap) {
                     const int residue_index = cached_residue.cmap_interaction.residue_index;
                     const unsigned int cmap_type_index = cached_residue.cmap_interaction.cmap_type_index;

                     const double phi = (*(this->chain))[residue_index].get_phi();
                     const double psi = (*(this->chain))[residue_index].get_psi();

                     energy_sum += charmm_cmap::cmap_energy(phi, psi, cmap_type_index, this->cmap_data);
               }
          }

          return energy_sum;

     }

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmmBondedCached(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm-bonded-cached", settings, random_number_engine) {

         this->none_move = false;
         setup_caches();
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmmBondedCached(const TermCharmmBondedCached &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {

          this->none_move = false;
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

          // By default don't treat energy evaluation as a none move
          this->none_move = false;

          if (move_info) {
 
               // This is a none move 
               if (move_info->modified_angles.empty() == true) {
 
                   // Notify accept/reject functions that this was a none_move
                   this->none_move = true;
 
                   // Return energy.
                   return this->energy_new * charmm_constants::KJ_TO_KCAL;
 
               // Not a none move
               } else {
 
                   // If these are set explicitly by the move, read these here.
                   this->start_index = std::max(0, move_info->modified_positions_start - 1);
                   this->end_index = std::min(move_info->modified_positions_end + 1, this->chain->size() - 1);
 
               }
          }

          // Local delta energy required for OpenMP -- can't just write to this->energy_new.
          double delta_energy_local = 0.0;

          // #pragma omp parallel for reduction(+:delta_energy_local) schedule(static)
          for (int i = this->start_index; i < this->end_index+1; i ++) {

               const double residue_energy = 
                    calculate_cached_residue_energy(this->bonded_cached_residues[i]);

               this->bonded_cached_residues[i].energy_new = residue_energy;

               delta_energy_local += this->bonded_cached_residues[i].energy_new
                                   - this->bonded_cached_residues[i].energy_old;
          }

          // Add energy delta
          this->energy_new  += delta_energy_local;

          // Return energy (and convert from kJ to kcal)
          return this->energy_new * charmm_constants::KJ_TO_KCAL;
     }


    //! Accept move and backup energies
     void accept() {

        if (this->none_move == false) {
            // If move is accepted, backup energies in all pairs that were recomputed
            for (int i = this->start_index; i < this->end_index+1; i ++) {

                this->bonded_cached_residues[i].energy_old
                    = this->bonded_cached_residues[i].energy_new;
            }
            this->energy_old = this->energy_new;
        }
    }


    //! Reject move and roll-back energies
    void reject() {

        if (this->none_move == false) {

            // If move is accepted, restore energies in all pairs that were recomputed
            for (int i = this->start_index; i < this->end_index+1; i ++) {

                this->bonded_cached_residues[i].energy_new
                    = this->bonded_cached_residues[i].energy_old;
            }
            this->energy_new = this->energy_old;
        }
    }


};

} // End namespace phaistos
#endif
