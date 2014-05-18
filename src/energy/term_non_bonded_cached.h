// term_vdw.h --- Van der Waals interaction energy term
// Copyright (C) 2014 Anders S. Christensen
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

#ifndef TERM_GROMACS_NON_BONDED_CACHED_H
#define TERM_GROMACS_NON_BONDED_CACHED_H

#include <string>

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "charmm22_parser.h"

namespace phaistos {

//! Gromacs van der Waals interaction term
class TermGromacsNonBondedCached: public EnergyTermCommon<TermGromacsNonBondedCached, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermGromacsNonBondedCached, ChainFB> EnergyTermCommon;

     //! Number of interactions calculated
     int counter;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     struct CachedResiduePair {

        double energy_old;
        double energy_new;

        std::vector<NonBondedPair> pairs;

     };


     std::vector<NonBondedParameter> non_bonded_parameters;
     std::vector<NonBonded14Parameter> non_bonded_14_parameters;
     std::vector<NonBondedPair> non_bonded_pairs;
     std::vector< std::vector<CachedResiduePair> > cached_residue_pairs;
     double total_energy;
     double total_energy_old;
     std::vector<std::vector<unsigned int> > cache_indexes;


     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermGromacsNonBondedCached(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "gromacs-non-bonded-cached", settings, random_number_engine) {

              std::string non_bonded_filename = "/home/andersx/phaistos_dev/modules/gromacs/src/energy/charmm22_cmap/charmm22_vdw.itp";
              non_bonded_parameters = read_nonbonded_parameters(non_bonded_filename);

              std::string non_bonded_14_filename = "/home/andersx/phaistos_dev/modules/gromacs/src/energy/charmm22_cmap/charmm22_vdw14.itp";
              non_bonded_14_parameters = read_nonbonded_14_parameters(non_bonded_14_filename);

              non_bonded_pairs = generate_non_bonded_pairs_cached(this->chain,
                                                           non_bonded_parameters,
                                                           non_bonded_14_parameters);

            for (int i = 0; i < this->chain->size(); i++) {
                std::vector<CachedResiduePair> temp_pair_vector;
                this->cached_residue_pairs.push_back(temp_pair_vector);
                for (int j = 0; j < this->chain->size(); j++) {
                    CachedResiduePair temp_pair;
                    temp_pair.energy_old = 0.0;
                    temp_pair.energy_new = 0.0;
                    this->cached_residue_pairs[i].push_back(temp_pair);
                }
            }

            for (unsigned int i = 0; i < non_bonded_pairs.size(); i++) {
                NonBondedPair pair = non_bonded_pairs[i];

                int residue1_index = (pair.atom1)->residue->index;
                int residue2_index = (pair.atom2)->residue->index;

                if (residue2_index > residue1_index) {
                    this->cached_residue_pairs[residue2_index][residue1_index].pairs.push_back(pair);
                } else {
                    this->cached_residue_pairs[residue1_index][residue2_index].pairs.push_back(pair);
                }
            }

            total_energy = 0.0;
            total_energy_old = 0.0;

            for (int i = 0; i < this->chain->size(); i++) {
                for (int j = 0; j < this->chain->size(); j++) {

                    this->cached_residue_pairs[i][j].energy_old = 0.0;
                    this->cached_residue_pairs[i][j].energy_new = 0.0;
                    // if (j > i) continue;

                    for (unsigned int k = 0; k < this->cached_residue_pairs[i][j].pairs.size(); k++) {

                        NonBondedPair pair = this->cached_residue_pairs[i][j].pairs[k];

                        const double r_sq = ((pair.atom1)->position - (pair.atom2)->position).norm_squared();
                        const double inv_r_sq = 100.0 / (r_sq); //shift to nanometers
                        const double inv_r_sq6 = inv_r_sq * inv_r_sq * inv_r_sq;
                        const double inv_r_sq12 = inv_r_sq6 * inv_r_sq6;
                        const double vdw_energy_temp = (pair.c12 * inv_r_sq12 - pair.c6 * inv_r_sq6) / 4.184;

                        const double coul_energy_temp = pair.qq * sqrt(inv_r_sq)/ 4.184;

                        this->cached_residue_pairs[i][j].energy_old += vdw_energy_temp + coul_energy_temp;
                        this->cached_residue_pairs[i][j].energy_new += vdw_energy_temp + coul_energy_temp;

                        total_energy     += vdw_energy_temp + coul_energy_temp;
                        total_energy_old += vdw_energy_temp + coul_energy_temp;

                    }
                }
            }

            std::cout << "Total constructor energy " << total_energy << std::endl;


            cache_indexes = get_cache_indexes(0, this->chain->size() -1);

     }


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


     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermGromacsNonBondedCached(const TermGromacsNonBondedCached &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            counter(other.counter),
            non_bonded_pairs(other.non_bonded_pairs),
            cached_residue_pairs(other.cached_residue_pairs),
            total_energy(other.total_energy),
            total_energy_old(other.total_energy_old),
            cache_indexes(other.cache_indexes) {
                std::cout << "Executed copy constructor" << std::endl;
            }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

         // Indexes of first and last residue for which the
         // position of atoms have changed since last move.
         unsigned int start_index = 0;
         unsigned int end_index = this->chain->size() - 1;

         // If these are set explicitly by the move, read these here.
         if (move_info) {
             start_index = move_info->modified_positions_start;
             end_index = move_info->modified_positions_end - 1;
         }

        // Get the indexes of all pairs of residues which must be recomputed
        this->cache_indexes = get_cache_indexes(start_index, end_index);

        // Local delta energy required for OpenMP.
        double delta_energy_local = 0.0;

        // Loop over all pairs which must be recomputed
        #pragma omp parallel for reduction(+:delta_energy_local) schedule(static)
        for (unsigned int k = 0; k < this->cache_indexes.size(); k++) {

            // Local index variables for residues i and j.
            unsigned int i = this->cache_indexes[k][0];
            unsigned int j = this->cache_indexes[k][1];

            // Reset interaction energy of the residue pair ij.
            this->cached_residue_pairs[i][j].energy_new = 0.0;

            // Loop over all pairs of atoms, k, in the residue pair ij.
            for (unsigned int k = 0; k < this->cached_residue_pairs[i][j].pairs.size(); k++) {

                // Make local copy of atom pair k.
                NonBondedPair pair = this->cached_residue_pairs[i][j].pairs[k];

                const double r_sq = ((pair.atom1)->position - (pair.atom2)->position).norm_squared();
                const double inv_r_sq = 100.0 / (r_sq); //shift to nanometers
                const double inv_r_sq6 = inv_r_sq * inv_r_sq * inv_r_sq;
                const double inv_r_sq12 = inv_r_sq6 * inv_r_sq6;

                const double vdw_energy_temp = (pair.c12 * inv_r_sq12 - pair.c6 * inv_r_sq6) / 4.184;
                const double coul_energy_temp = pair.qq * sqrt(inv_r_sq)/ 4.184;

                // Add diatomic contribution to interaction energy of the residue pair ij.
                this->cached_residue_pairs[i][j].energy_new += vdw_energy_temp + coul_energy_temp;
            }

            // Compute delta energy for the residue pair ij (I.e. subtract old, add new)
            delta_energy_local += this->cached_residue_pairs[i][j].energy_new
                                - this->cached_residue_pairs[i][j].energy_old;

        }

        // Add delta energy for the move to total energy.
        total_energy += delta_energy_local;

        // Return energy.
        return total_energy;

     }



    void accept() {
        for (unsigned int k = 0; k < this->cache_indexes.size(); k++) {
            unsigned int i = this->cache_indexes[k][0];
            unsigned int j = this->cache_indexes[k][1];
            this->cached_residue_pairs[i][j].energy_old = this->cached_residue_pairs[i][j].energy_new;
        }
        total_energy_old = total_energy;
    }

    void reject() {
        for (unsigned int k = 0; k < this->cache_indexes.size(); k++) {
            unsigned int i = this->cache_indexes[k][0];
            unsigned int j = this->cache_indexes[k][1];
            this->cached_residue_pairs[i][j].energy_new = this->cached_residue_pairs[i][j].energy_old;
        }
        total_energy = total_energy_old;
    }

};

}

#endif
