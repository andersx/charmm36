// non_bonded_cached.h --- Cached CHARMM36/EEF1-SB Coulomb, van der Waals and implicit solvent collective term
// Copyright (C) 2014 Sandro Bottaro, Anders S. Christensen
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

#ifndef TERM_CHARMM_NON_BONDED_CACHED_H
#define TERM_CHARMM_NON_BONDED_CACHED_H

#include <string>

#include <boost/type_traits/is_base_of.hpp>
#include <boost/tokenizer.hpp>

#include "energy/energy_term.h"

#include "parsers/topology_parser.h"
#include "parsers/eef1_sb_parser.h"
#include "constants.h"
#include "parameters/vdw14_itp.h"
#include "parameters/vdw_itp.h"
#include "parameters/solvpar_17_inp.h"

namespace phaistos {

class TermCharmmNonBondedCached: public EnergyTermCommon<TermCharmmNonBondedCached, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmmNonBondedCached, ChainFB> EnergyTermCommon;

     //! Struct that holds information that need to be computed if
     //! interactions between two residues change.
     struct CachedResidueInteraction {

        //! Energy before latest move
        double energy_old;

        //! Energy after latest move
        double energy_new;

        //! A list of all interactions that need to be computed
        std::vector<topology::NonBondedInteraction> interactions;

     };

     //! Lookup table containing interactions that need to be recomputed.
     //! The size is L x L where L is the length of the chain. Only the upper triangle
     //! is filled.
     std::vector< std::vector<CachedResidueInteraction> > cached_residue_interactions;

     //! The total energy after the latest move
     double total_energy;

     //! The total energy before the latest move
     double total_energy_old;

     //! List of indexes in the cache that need to be recomputed after last move
     std::vector<std::vector<unsigned int> > cache_indexes;

     bool none_move;

     double dGref_total;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;


     //! Calculate the energy of a diatomic interaction
     double calculate_interaction_energy(topology::NonBondedInteraction &interaction) {

          const double r_sq = ((interaction.atom1)->position - (interaction.atom2)->position).norm_squared();

          const double inv_r_sq = 1.0 / r_sq;
          const double inv_r_sq6 = inv_r_sq * inv_r_sq * inv_r_sq * charmm_constants::NM6_TO_ANGS6; //shift to nanometers^6
          const double inv_r_sq12 = inv_r_sq6 * inv_r_sq6;

          const double vdw_energy = (interaction.c12 * inv_r_sq12 - interaction.c6 * inv_r_sq6);
          const double coul_energy = interaction.qq * inv_r_sq * charmm_constants::TEN_OVER_ONE_POINT_FIVE;

          double eef1_sb_energy = 0.0;

          // If the pair has a contribution to EEF1-SB
          if ((interaction.do_eef1) && (r_sq < 81.0)) {

              // From Sandro's code
              const double r_ij = std::sqrt(r_sq);

              double R_min_i = interaction.R_vdw_1;
              double R_min_j = interaction.R_vdw_2;

              double lambda_i = interaction.lambda1;
              double lambda_j = interaction.lambda2;

              const double arg_ij = std::fabs((r_ij - R_min_i)/lambda_i);
              const double arg_ji = std::fabs((r_ij - R_min_j)/lambda_j);

              int bin_ij = int(arg_ij*100);
              int bin_ji = int(arg_ji*100);

              double exp_ij = 0.0;
              double exp_ji = 0.0;

              if (bin_ij < 350) exp_ij = charmm_constants::EXP_EEF1[bin_ij];
              if (bin_ji < 350) exp_ji = charmm_constants::EXP_EEF1[bin_ji];

              double cont_ij = -interaction.fac_12*exp_ij * inv_r_sq;
              double cont_ji = -interaction.fac_21*exp_ji * inv_r_sq;

              // Add to local EEF1-SB energy
              eef1_sb_energy += cont_ij + cont_ji;

          }

          return (vdw_energy + coul_energy) * charmm_constants::KJ_TO_KCAL + eef1_sb_energy;
     }


     //! Code that needs to be setup explicitly in the constructor, and also
     //! explicitly copy-constructor. Sets up the cache and initial energies.
     void setup_caches() {


          std::vector<double> dGref;
          std::vector< std::vector<double> > factors;
          std::vector<double> vdw_radii;
          std::vector<double> lambda;
          std::map<std::string, unsigned int> eef1_atom_type_index_map;

          initialize(dGref, factors, vdw_radii, lambda, eef1_atom_type_index_map);


          std::vector<topology::NonBondedParameter> non_bonded_parameters
              = topology::read_nonbonded_parameters(charmm_constants::vdw_itp);

          std::vector<topology::NonBonded14Parameter> non_bonded_14_parameters
              = topology::read_nonbonded_14_parameters(charmm_constants::vdw14_itp);

          std::vector<topology::NonBondedInteraction> non_bonded_interactions
              = topology::generate_non_bonded_interactions_cached(this->chain,
                                                                  non_bonded_parameters,
                                                                  non_bonded_14_parameters,
                                                                  dGref,
                                                                  factors,
                                                                  vdw_radii,
                                                                  lambda,
                                                                  eef1_atom_type_index_map);

            std::cout << non_bonded_interactions.size() << std::endl;
            this->dGref_total = 0.0;

            for (AtomIterator<ChainFB, definitions::ALL> it(*this->chain); !it.end(); ++it) {

                Atom *atom = &*it;
                std::string atom_type = eef1_sb_parser::get_atom_type(atom);

                unsigned int index = eef1_atom_type_index_map[atom_type];

                this->dGref_total += dGref[index];
            }

            // Fill up cached matrix and set energy to zero
            for (int i = 0; i < this->chain->size(); i++) {

                // Dummy residue pair
                std::vector<CachedResidueInteraction> temp_interaction_vector;
                this->cached_residue_interactions.push_back(temp_interaction_vector);

                for (int j = 0; j < this->chain->size(); j++) {

                    // Initialize residue pair and push back
                    CachedResidueInteraction temp_interaction;
                    temp_interaction.energy_old = 0.0;
                    temp_interaction.energy_new = 0.0;
                    this->cached_residue_interactions[i].push_back(temp_interaction);
                }
            }

            // Fill atom pairs in cache matrix
            for (unsigned int i = 0; i < non_bonded_interactions.size(); i++) {

                topology::NonBondedInteraction interaction = non_bonded_interactions[i];

                // Get the residue indexes of the atoms involved in this pair
                int residue1_index = (interaction.atom1)->residue->index;
                int residue2_index = (interaction.atom2)->residue->index;

                // Make sure we only fill out the upper triangle of the cache matrix
                if (residue2_index > residue1_index) {
                    this->cached_residue_interactions[residue2_index][residue1_index].interactions.push_back(interaction);
                } else {
                    this->cached_residue_interactions[residue1_index][residue2_index].interactions.push_back(interaction);
                }
            }

            // Initialize total energies
            this->total_energy = this->dGref_total;
            this->total_energy_old = this->dGref_total;

            // Calculate energy for each cache matrix element
            for (int i = 0; i < this->chain->size(); i++) {
                for (int j = 0; j < this->chain->size(); j++) {

                    // Reset matrix element energies
                    this->cached_residue_interactions[i][j].energy_old = 0.0;
                    this->cached_residue_interactions[i][j].energy_new = 0.0;
                    // if (j > i) continue;

                    // Sum over all interactions in that matrix element
                    for (unsigned int k = 0; k < this->cached_residue_interactions[i][j].interactions.size(); k++) {

                        topology::NonBondedInteraction interaction = this->cached_residue_interactions[i][j].interactions[k];

                        double interaction_energy = calculate_interaction_energy(interaction);

                        // Add to matrix element
                        this->cached_residue_interactions[i][j].energy_old += interaction_energy;
                        this->cached_residue_interactions[i][j].energy_new += interaction_energy;

                        //Add to toal energies
                        this->total_energy     += interaction_energy;
                        this->total_energy_old += interaction_energy;
                    }
                }
            }

            std::cout << "Total constructor energy " << this->total_energy << std::endl;


            this->cache_indexes = get_cache_indexes(0, this->chain->size() -1);
     }


     // Returns the pairs indexes for the residues which need recomputation after
     // an MC move has been carried out.
     std::vector<std::vector<unsigned int> > get_cache_indexes(const unsigned int start,
                                                               const unsigned int end) {
        using namespace vector_utils;
        std::vector<std::vector<unsigned int> > indexes;

        unsigned int last = (unsigned int)this->chain->size();

        for (unsigned int i = start; i < end + 1; i++) {
            for (unsigned int j = 0; j < last; j++) {

                if (j < start) {
                    indexes.push_back(make_vector(i, j));
                } else if (j > end) {
                    indexes.push_back(make_vector(j, i));
                } else if (i >= j) {
                    indexes.push_back(make_vector(i, j));
                }
            }
        }

        return indexes;
     }


     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmmNonBondedCached(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm-non-bonded-cached", settings, random_number_engine) {

          this->none_move = false;
          setup_caches();
     }


     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmmNonBondedCached(const TermCharmmNonBondedCached &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {

          this->none_move = false;
          setup_caches();
     }


     // Big long initialize code from Wouter/Sandro
     void initialize(std::vector<double> &dGref,
                     std::vector< std::vector<double> > &factors,
                     std::vector<double> &vdw_radii,
                     std::vector<double> &lambda,
                     std::map<std::string, unsigned int> &eef1_atom_type_index_map) {

          // Read parameter file
          std::vector<std::string> atoms;
          std::vector< std::vector<double> > params;
          std::string elem_str;

          //! Useful constant
          const double two_pi_3_2 = 2.0*M_PI*sqrt(M_PI);
          const double phys_t = 298.15;

          std::istringstream f_h(charmm_constants::solvpar_inp);

          std::string line;
          std::vector<double> pp;

          unsigned int eef1_atom_type_index = 0;

          while (getline(f_h,line)) {

               //Skip comments
               if (line.find("!") == 0)
                    continue;

               int index = 0;

               boost::char_separator<char> sep(" ");
               boost::tokenizer<boost::char_separator<char> > tok(line,sep);
               for(boost::tokenizer<boost::char_separator<char> >::iterator beg=tok.begin(); beg!=tok.end();++beg){
                    index++;
                    if(index==1){

                         atoms.push_back((*beg).c_str());

                         eef1_atom_type_index_map[(*beg).c_str()] = eef1_atom_type_index;
                         eef1_atom_type_index++;

                    } else {
                         double elem_val=strtod( (*beg).c_str(),NULL );
                         pp.push_back(elem_val);
                    }
               }
               params.push_back(pp);
               pp.clear();
          }

          //Create pairwise paramters in a atomtype*atomtype lookup table
          //double t = settings.temp;
          double t = 298.15;
          double dt = t-phys_t;
          for (unsigned int i=0; i< atoms.size(); i++){
               //dGref(t) = dGref_t0 - (dH_t0-dGref_t0)*(dt/t_0) - dCp(t*log(t/t_0) -dt)
               double cont_1 = (params[i][1]);
               double cont_2 = (params[i][3] - (params[i][1]))*(dt/phys_t);
               double cont_3 = 0.0;
               if(std::fabs(dt)>0.00001){
                    double a = std::log(t/phys_t);
                      cont_3 = (t*a - dt);
               }
               double dGref_i =  cont_1 - cont_2 - cont_3;
               dGref.push_back(dGref_i);
               vdw_radii.push_back(params[i][6]);
               lambda.push_back(params[i][5]);
          }

          for (unsigned int i=0; i<atoms.size(); i++){

               // dGfree_i(t) = (dGref(t)/dGref_t0)*dGfree_t0
               double dGfree_i = 0.0;
               if(std::fabs(params[i][1])>0.00001){
                    dGfree_i = (dGref[i]/(params[i][1]))*(params[i][2]);
               }

               std::vector<double> factors_i;
               for (unsigned int j=0; j<atoms.size(); j++){

                    //Factor = dGfree_i*V_j/(2*pi*sqrt(pi)*lambda_i)
                    double factor = 0.0;

                    if(std::fabs(params[i][5])>0.00001){
                         factor = (dGfree_i*params[j][0])/(two_pi_3_2*params[i][5]);
                    }

                    factors_i.push_back(factor);
               }
               factors.push_back(factors_i);
          }

     }


     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

         // Indexes of first and last residue for which the
         // position of atoms have changed since last move.
         unsigned int start_index = 0;
         unsigned int end_index = this->chain->size() - 1;

         this->none_move = false;

         if (move_info) {

            // This is a none move
            if (move_info->modified_angles.empty() == true) {

                  // Notify accept/reject functions that this was a none_move
                  this->none_move = true;

                  // Return energy.
                  return this->total_energy;

             // Not a none move
             } else {

                // If these are set explicitly by the move, read these here.
                start_index = move_info->modified_positions_start;
                end_index = move_info->modified_positions_end - 1;

            }
        }

        // Get the indexes of all pairs of residues which must be recomputed
        this->cache_indexes = get_cache_indexes(start_index, end_index);

        // Local delta energy required for OpenMP.
        double delta_energy_local = 0.0;

        // Loop over all pairs which must be recomputed
        // #pragma omp parallel for reduction(+:delta_energy_local) schedule(static)
        for (unsigned int k = 0; k < this->cache_indexes.size(); k++) {

            // Local index variables for residues i and j.
            unsigned int i = this->cache_indexes[k][0];
            unsigned int j = this->cache_indexes[k][1];

            // Reset interaction energy of the residue pair ij.
            this->cached_residue_interactions[i][j].energy_new = 0.0;

            // Loop over all pairs of atoms, k, in the residue pair ij.
            for (unsigned int k = 0; k < this->cached_residue_interactions[i][j].interactions.size(); k++) {

                // This is the loop where the majority of the time is spent. Feel free to optimize!

                // Make local copy of atom pair k.
                topology::NonBondedInteraction interaction = this->cached_residue_interactions[i][j].interactions[k];

                const double r2 = ((interaction.atom1)->position - (interaction.atom2)->position).norm_squared();

                const double inv_r2 = 1.0 / r2; // convert to nanometers
                const double inv_r6 = inv_r2 * inv_r2 * inv_r2 * charmm_constants::NM6_TO_ANGS6;

                // Add vdw and coulomb energy (using nm and kJ).
                this->cached_residue_interactions[i][j].energy_new += (interaction.c12 * inv_r6 - interaction.c6) * inv_r6                    // VDW energy
                                                                       + interaction.qq * inv_r2 * charmm_constants::TEN_OVER_ONE_POINT_FIVE;   // Coulomb energy

                // If the pair has a contribution to EEF1-SB solvation term
                if ((interaction.do_eef1) && (r2 < 81.0)) {

                    // From Sandro's code -- this bit is in angstrom and kcal.
                    const double r_ij = std::sqrt(r2);

                    const double arg_ij = std::fabs((r_ij - interaction.R_vdw_1)/interaction.lambda1);
                    const double arg_ji = std::fabs((r_ij - interaction.R_vdw_2)/interaction.lambda2);

                    const int bin_ij = int(arg_ij*100);
                    const int bin_ji = int(arg_ji*100);

                    double exp_ij = 0.0;
                    double exp_ji = 0.0;

                    if (bin_ij < 350) exp_ij = charmm_constants::EXP_EEF1[bin_ij];
                    if (bin_ji < 350) exp_ji = charmm_constants::EXP_EEF1[bin_ji];

                    // Add solvation energy (in kcal, so convert to kJ)
                    this->cached_residue_interactions[i][j].energy_new -= (interaction.fac_12*exp_ij + interaction.fac_21*exp_ji) * inv_r2 * charmm_constants::KCAL_TO_KJ;

                }

            }

            // Energies are summed in kJ, so convert to kcal now.
            this->cached_residue_interactions[i][j].energy_new *= charmm_constants::KJ_TO_KCAL;

            // Compute delta energy for the residue pair ij (I.e. subtract old, add new)
            delta_energy_local += this->cached_residue_interactions[i][j].energy_new
                                - this->cached_residue_interactions[i][j].energy_old;

        }

        // Select update scheme for total energy (necessary for precision when the change in energy is large)
        // If the energy-difference for the move is very large, the energy update might break double precision.
        // Required accuracy is around 1e-7, so if dE is > 1e6, then do a full summation to be safe.
        if (std::fabs(delta_energy_local) > 1e6) {

            // Get reference solvation energy of system
            this->total_energy = this->dGref_total;

            // Get the indexes of all pairs of residues which must be recomputed
            std::vector<std::vector<unsigned int> > all_indexes = get_cache_indexes(0, this->chain->size() - 1);

            for (unsigned int k = 0; k < all_indexes.size(); k++) {

                // Local index variables for residues i and j.
                unsigned int i = all_indexes[k][0];
                unsigned int j = all_indexes[k][1];
                this->total_energy += this->cached_residue_interactions[i][j].energy_new;
            }

        // If the energy difference is small, then just add the delta energy.
        } else {

            // Add delta energy for the move to total energy.
            this->total_energy += delta_energy_local;
        }

        // Return energy.
        return this->total_energy;

    }


    //! Accept move and backup energies
    void accept() {

        if (this->none_move == false) {

            // If move is accepted, backup energies in all pairs that were recomputed
            for (unsigned int k = 0; k < this->cache_indexes.size(); k++) {
                unsigned int i = this->cache_indexes[k][0];
                unsigned int j = this->cache_indexes[k][1];
                this->cached_residue_interactions[i][j].energy_old = this->cached_residue_interactions[i][j].energy_new;
            }

            //Backup total energy
            this->total_energy_old = this->total_energy;
        }
    }


    //! Reject move and roll-back energies
    void reject() {

        if (this->none_move == false) {

            // If move is accepted, restore energies in all pairs that were recomputed
            for (unsigned int k = 0; k < this->cache_indexes.size(); k++) {
                unsigned int i = this->cache_indexes[k][0];
                unsigned int j = this->cache_indexes[k][1];
                this->cached_residue_interactions[i][j].energy_new = this->cached_residue_interactions[i][j].energy_old;
            }

            // Restore total energy
            this->total_energy = this->total_energy_old;
        }
    }


};

} // End namespace phaistos

#endif
