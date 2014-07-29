// non_bonded.h --- Coulomb, vdw and excl vol solvent collective term
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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Phaistos. If not, see <http://www.gnu.org/licenses/>.
//

#ifndef TERM_CHARMM36_NON_BONDED_H
#define TERM_CHARMM36_NON_BONDED_H

#include <string>

#include <boost/tokenizer.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "parsers/topology_parser.h"
#include "parsers/eef1_sb_parser.h"
#include "constants.h"
#include "parameters/vdw14_itp.h"
#include "parameters/vdw_itp.h"
#include "parameters/solvpar_17_inp.h"

namespace phaistos {

class TermCharmm36NonBonded: public EnergyTermCommon<TermCharmm36NonBonded, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36NonBonded, ChainFB> EnergyTermCommon;

     std::vector<topology::NonBondedInteraction> non_bonded_interactions;
     double dGref_total;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;


     void setup_interactions() {


          std::vector<double> dGref;
          std::vector< std::vector<double> > factors;
          std::vector<double> vdw_radii;
          std::vector<double> lambda;
          std::map<std::string, unsigned int> eef1_atom_type_index_map;

          initialize(dGref, factors, vdw_radii, lambda, eef1_atom_type_index_map);


          std::vector<topology::NonBondedParameter> non_bonded_parameters
              = topology::read_nonbonded_parameters(charmm36_constants::vdw_itp);

          std::vector<topology::NonBonded14Parameter> non_bonded_14_parameters
              = topology::read_nonbonded_14_parameters(charmm36_constants::vdw14_itp);

          this->non_bonded_interactions = topology::generate_non_bonded_interactions_cached(this->chain,
                                                                  non_bonded_parameters,
                                                                  non_bonded_14_parameters,
                                                                  dGref,
                                                                  factors,
                                                                  vdw_radii,
                                                                  lambda,
                                                                  eef1_atom_type_index_map);

            std::cout << this->non_bonded_interactions.size() << std::endl;
            this->dGref_total = 0.0;

            for (AtomIterator<ChainFB, definitions::ALL> it(*this->chain); !it.end(); ++it) {

                Atom *atom = &*it;
                std::string atom_type = eef1_sb_parser::get_atom_type(atom);

                unsigned int index = eef1_atom_type_index_map[atom_type];

                this->dGref_total += dGref[index];
            }
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

          std::istringstream f_h(charmm36_constants::solvpar_inp);

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

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36NonBonded(ChainFB *chain,
                    const Settings &settings = Settings(),
                    RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-non-bonded", settings, random_number_engine) {

        setup_interactions();
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36NonBonded(const TermCharmm36NonBonded &other,
                 RandomNumberEngine *random_number_engine,
                 int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {

        setup_interactions();

        }


     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {


          double energy_sum = this->dGref_total * charmm36_constants::KCAL_TO_KJ;

          for (unsigned int i = 0; i < this->non_bonded_interactions.size(); i++) {

                // This is the loop where the majority of the time is spent. Feel free to optimize!

                // Make local copy of atom pair k.
                topology::NonBondedInteraction interaction = this->non_bonded_interactions[i];

                const double r2 = ((interaction.atom1)->position - (interaction.atom2)->position).norm_squared();

                const double inv_r2 = 1.0 / r2; // convert to nanometers
                const double inv_r6 = inv_r2 * inv_r2 * inv_r2 * charmm36_constants::NM6_TO_ANGS6;

                // Add vdw and coulomb energy (using nm and kJ).
                energy_sum += (interaction.c12 * inv_r6 - interaction.c6) * inv_r6     // VDW energy
                    + interaction.qq * inv_r2 * charmm36_constants::TEN_OVER_ONE_POINT_FIVE;   // Coulomb energy

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

                    if (bin_ij < 350) exp_ij = charmm36_constants::EXP_EEF1[bin_ij];
                    if (bin_ji < 350) exp_ji = charmm36_constants::EXP_EEF1[bin_ji];

                    // Add solvation energy (in kcal, so convert to kJ)
                    energy_sum -= (interaction.fac_12*exp_ij + interaction.fac_21*exp_ji) * inv_r2 * charmm36_constants::KCAL_TO_KJ;

                }


          }

          return energy_sum * charmm36_constants::KJ_TO_KCAL;
     }
};

}

#endif
