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

     struct CachedResiduePair {

        double energy_old;
        double energy_new;

        std::vector<NonBondedPair> pairs;

     };

     double total_energy;
     double total_energy_old;

     std::vector<NonBondedParameter> non_bonded_parameters;
     std::vector<NonBonded14Parameter> non_bonded_14_parameters;
     std::vector<NonBondedPair> non_bonded_pairs;

     std::vector<std::vector<unsigned int> > cache_indexes;
     std::vector< std::vector<CachedResiduePair> > cached_residue_pairs;


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


            for (unsigned int i = 0; i < this->chain->size(); i++) {
                std::vector<CachedResiduePair> temp_pair_vector;
                cached_residue_pairs.push_back(temp_pair_vector);
                for (unsigned int j = 0; j < this->chain->size(); j++) {

                    CachedResiduePair temp_pair;
                    cached_residue_pairs[i].push_back(temp_pair);

                }
            }
            for (unsigned int i = 0; i < non_bonded_pairs.size(); i++) {
                NonBondedPair pair = non_bonded_pairs[i];

                int residue1_index = (pair.atom1)->residue->index;
                int residue2_index = (pair.atom2)->residue->index;

                if (residue2_index > residue1_index) {
                    cached_residue_pairs[residue2_index][residue1_index].pairs.push_back(pair);
                }
                else {
                    cached_residue_pairs[residue1_index][residue2_index].pairs.push_back(pair);

                }

            }



            for (unsigned int i = 0; i < this->chain->size(); i++) {
                for (unsigned int j = 0; j < this->chain->size(); j++) {

                    cached_residue_pairs[i][j].energy_old = 0.0;
                    cached_residue_pairs[i][j].energy_new = 0.0;
                    // if (j > i) continue;

                    for (unsigned int k = 0; k < cached_residue_pairs[i][j].pairs.size(); k++) {

                        NonBondedPair pair = cached_residue_pairs[i][j].pairs[k];

                        const double r_sq = ((pair.atom1)->position - (pair.atom2)->position).norm_squared();
                        const double inv_r_sq = 100.0 / (r_sq); //shift to nanometers
                        const double inv_r_sq6 = inv_r_sq * inv_r_sq * inv_r_sq;
                        const double inv_r_sq12 = inv_r_sq6 * inv_r_sq6;
                        const double vdw_energy_temp = (pair.c12 * inv_r_sq12 - pair.c6 * inv_r_sq6) / 4.184;

                        const double coul_energy_temp = pair.qq * sqrt(inv_r_sq)/ 4.184;

                        cached_residue_pairs[i][j].energy_old += vdw_energy_temp + coul_energy_temp;
                        cached_residue_pairs[i][j].energy_new += vdw_energy_temp + coul_energy_temp;

                        total_energy     += vdw_energy_temp + coul_energy_temp;
                        total_energy_old += vdw_energy_temp + coul_energy_temp;

                    }
                }
            }

            std::cout << "Total constructor energy " << total_energy << std::endl;

        cache_indexes = get_cache_indexes(0, this->chain->size() -1);

     }


     std::vector<std::vector<unsigned int> > get_cache_indexes(unsigned int start, unsigned int end) {

        using namespace vector_utils;
        std::vector<std::vector<unsigned int> > indexes;

        for (unsigned int i = start; i < end + 1; i++) {
            for (unsigned int j = 0; j < this->chain->size(); j++) {

                if (j < start) {
                    indexes.push_back(make_vector(i, j));
                } else if (j > end) {
                    indexes.push_back(make_vector(i, j));
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
     TermGromacsVdw(const TermGromacsVdw &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            counter(other.counter),
            non_bonded_pairs(other.non_bonded_pairs) {}

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

            for (unsigned int i = 0; i < this->chain->size(); i++) {
                for (unsigned int j = 0; j < this->chain->size(); j++) {

                    std::cout << "SIZE: " << i << "  " << j << "  " << this->cached_residue_pairs[i][j].pairs.size() << "  " << this->cached_residue_pairs[i][j].energy_old << std::endl;

                }
            }
         unsigned int start_index = 0;
         unsigned int end_index = this->chain->size() - 1;

         if (move_info) {
             start_index = move_info->modified_positions_start;
             end_index = move_info->modified_positions_end-1;
         }


        cache_indexes = get_cache_indexes(start_index, end_index);

        total_energy = total_energy_old;

        for (unsigned int k = 0; k < cache_indexes.size(); k++) {

            unsigned int i = cache_indexes[k][0];
            unsigned int j = cache_indexes[k][1];

            printf("Indexes:  i= %4d  j=%4d\n", i, j);
        }
        for (unsigned int k = 0; k < cache_indexes.size(); k++) {

            unsigned int i = cache_indexes[k][0];
            unsigned int j = cache_indexes[k][1];

            printf("Indexes:  i= %4d  j=%4d\n", i, j);

            // cached_residue_pairs[i][j].energy_new = 0.0;
            std::cout << cached_residue_pairs.size() << std::endl; //= 0.0;
            std::cout << cached_residue_pairs[i].size() << std::endl; //= 0.0;
            // std::cout << cached_residue_pairs[i][j].energy_new << std::endl; //= 0.0;

            printf("Did update energy.\n");

            for (unsigned int k = 0; k < cached_residue_pairs[i][j].pairs.size(); k++) {

                printf("Trying:  k= %4d \n", k);

                NonBondedPair pair = cached_residue_pairs[i][j].pairs[k];

                const double r_sq = ((pair.atom1)->position - (pair.atom2)->position).norm_squared();
                const double inv_r_sq = 100.0 / (r_sq); //shift to nanometers
                const double inv_r_sq6 = inv_r_sq * inv_r_sq * inv_r_sq;
                const double inv_r_sq12 = inv_r_sq6 * inv_r_sq6;
                const double vdw_energy_temp = (pair.c12 * inv_r_sq12 - pair.c6 * inv_r_sq6) / 4.184;

                const double coul_energy_temp = pair.qq * sqrt(inv_r_sq)/ 4.184;

                cached_residue_pairs[i][j].energy_new += vdw_energy_temp + coul_energy_temp;

            }

            total_energy += cached_residue_pairs[i][j].energy_new - cached_residue_pairs[i][j].energy_old;

        }

          double energy_sum = 0.0;

          double vdw_energy = 0.0;
          double coul_energy = 0.0;


          double vdw14_energy = 0.0;
          double coul14_energy = 0.0;
          std::cout << "herelater" << std::endl;

          #pragma omp parallel for reduction(+:energy_sum) schedule(static)
          for (unsigned int i = 0; i < non_bonded_pairs.size(); i++) {

               NonBondedPair pair = non_bonded_pairs[i];

               const double r_sq = ((pair.atom1)->position - (pair.atom2)->position).norm_squared();
               const double inv_r_sq = 100.0 / (r_sq); //shift to nanometers
               const double inv_r_sq6 = inv_r_sq * inv_r_sq * inv_r_sq;
               const double inv_r_sq12 = inv_r_sq6 * inv_r_sq6;
               const double vdw_energy_temp = (pair.c12 * inv_r_sq12 - pair.c6 * inv_r_sq6) / 4.184;

               const double coul_energy_temp = pair.qq * sqrt(inv_r_sq)/ 4.184;

               energy_sum += vdw_energy_temp + coul_energy_temp;

          }

          std::cout << "uncached energy sum " << energy_sum << std::endl;
          std::cout << "cached energy sum " << total_energy << std::endl;
          return energy_sum;
     }

    void accept() {
        // std::vector<std::vector<unsigned int> > cache_indexes = get_cache_indexes(start_index, end_index);

        for (unsigned int k = 0; k < cache_indexes.size(); k++) {

            unsigned int i = cache_indexes[k][0];
            unsigned int j = cache_indexes[k][1];

            cached_residue_pairs[i][j].energy_old = cached_residue_pairs[i][j].energy_new;


        }
        total_energy_old = total_energy;
    }

    void reject() {
        // std::vector<std::vector<unsigned int> > cache_indexes = get_cache_indexes(start_index, end_index);

        for (unsigned int k = 0; k < cache_indexes.size(); k++) {

            unsigned int i = cache_indexes[k][0];
            unsigned int j = cache_indexes[k][1];

            cached_residue_pairs[i][j].energy_new = cached_residue_pairs[i][j].energy_old;


        }
        total_energy = total_energy_old;

    }



};

}

#endif
