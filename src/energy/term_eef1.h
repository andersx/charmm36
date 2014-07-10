// term_eef1.h -- EEF1-SB solvation term.
// Copyright (C) 2009-2014 Sandro Bottaro, Anders S. Christensen
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

#ifndef TERM_CHARMM36_EEF1_H
#define TERM_CHARMM36_EEF1_H

#include <boost/type_traits/is_base_of.hpp>
#include <boost/tokenizer.hpp>
#include "energy/energy_term.h"
#include "parsers/eef1_sb_parser.h"

namespace phaistos {

//! partial eef1 interaction term
class TermCharmm36Eef1: public EnergyTermCommon<TermCharmm36Eef1, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36Eef1, ChainFB> EnergyTermCommon;

     //! Number of interactions in the last evaluation
     int counter;

protected:

     //! Lookup tables containing parameters
     std::vector<double> dGref;
     std::vector< std::vector<double> > factors;
     std::vector<double> vdw_radii;
     std::vector<double> lambda;

    std::map<std::string, unsigned int> eef1_atom_type_index_map;

public:


     //! Local settings class
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:

          //! Path to file containing solvation parameters
          std::string solvation_filename;

          //! Temperature
          double temp;

          //! Constructor
          Settings(std::string solvation_filename 
                  = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/parameters/solvpar_17.inp",double temp=298.15)
               : solvation_filename(solvation_filename),temp(temp) {}

          //! Output operator
          friend std::ostream &operator<<(std::ostream &o, const Settings &settings) {
               o << "solvation-filename:" << settings.solvation_filename << "\n";
               o << "temperature:" << settings.temp << "\n";
               o << static_cast<const EnergyTerm<ChainFB>::Settings>(settings);
               return o;
          }
     } settings;    //!< Local settings object

     //! Constructor
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36Eef1(ChainFB *chain,
              const Settings &settings=Settings(),
              RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-eef1", settings, random_number_engine) {

          initialize();

     }


     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36Eef1(const TermCharmm36Eef1 &other,
                     RandomNumberEngine *random_number_engine,
                     int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            counter(other.counter){
          initialize();
     }


     void initialize() {

          // Read parameter file
          std::vector<std::string> atoms;
          std::vector< std::vector<double> > params;
          std::string elem_str;

          //! Useful constant
          const double two_pi_3_2 = 2.0*M_PI*sqrt(M_PI);
          const double phys_t = 298.15;

          std::ifstream f_h((settings.solvation_filename).c_str());
          if (!f_h) {
               printf("# ERROR: EEF1 input: unable to open file '%s'\n", (settings.solvation_filename).c_str() );
               exit(1); // terminate with error
          }

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
          f_h.close();

          //Create pairwise paramters in a atomtype*atomtype lookup table
          double t = settings.temp;
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

     //! Evaluate eef1 interaction between 2 atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param chg1 Eef1 of atom1
     //! \param chg2 Eef1 of atom2
     //! \return Eef1 energy for atom pair
     double calculate_contribution(Atom *atom1, Atom *atom2, int index1, int index2) {

             return calc_eef1_energy(atom1,atom2,index1,index2);
     }

     //! Evaluate a eef1 interaction between two atoms
     //! \param atom1 First atom
     //! \param atom2 Second atom
     //! \param chg1 Eef1 of atom1
     //! \param chg2 Eef1 of atom2
     //! \return Eef1 energy for atom pair
     double calc_eef1_energy(Atom *atom1, Atom *atom2, int index1, int index2) {

          counter++;

          double r_ij = (atom1->position - atom2->position).norm();
          double r_ij_sq = r_ij*r_ij;
          double R_min_i = vdw_radii[index1];
          double R_min_j = vdw_radii[index2];

          double lambda_i = lambda[index1];
          double lambda_j = lambda[index2];

          double arg_ij = std::fabs((r_ij - R_min_i)/lambda_i);
          double arg_ji = std::fabs((r_ij - R_min_j)/lambda_j);

          // In CHARMM the exponential is not calculated explicitly
          // A lookup table is used instead
          int bin_ij = int(arg_ij*100);
          int bin_ji = int(arg_ji*100);

          double exp_ij = 0.0;
          double exp_ji = 0.0;

          if (bin_ij < 350) exp_ij = charmm36_constants::EXP_EEF1[bin_ij];
          if (bin_ji < 350) exp_ji = charmm36_constants::EXP_EEF1[bin_ji];

          double cont_ij = -factors[index1][index2]*exp_ij/r_ij_sq;
          double cont_ji = -factors[index2][index1]*exp_ji/r_ij_sq;

          // The code below is using the full std::exp function rather than lookup-table
          // double cont_ij = -factors[index1][index2]*std::exp(-(arg_ij*arg_ij))/r_ij_sq;
          // double cont_ji = -factors[index2][index1]*std::exp(-(arg_ji*arg_ji))/r_ij_sq;

          return (cont_ij+cont_ji);
     }

     //! Return index in parameter table corresponding to the atomtype
     int get_index(Atom *atom){

          std::string atom_type = eef1_sb_parser::get_atom_type(atom);
          unsigned int eef1_atom_type_index = eef1_atom_type_index_map[atom_type];

          return (int)eef1_atom_type_index;
     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energy_sum=0.0;
          counter=0;

          unsigned int i_atom = 0;

          double dGref_total = 0.0;

          // Iterate all the atom pairs on the chain
          for (AtomIterator<ChainFB, definitions::ALL> it1(*this->chain); !it1.end(); ++it1) {
               i_atom++;
               Atom *atom1 = &*it1;
               int index1 = get_index(atom1);

               if (atom1->mass == definitions::atom_h_weight) continue;

               energy_sum += dGref[index1];
               dGref_total += dGref[index1];

               unsigned int j_atom = 0;

               for (AtomIterator<ChainFB, definitions::ALL> it2(*this->chain); !it2.end(); ++it2) {
                    j_atom++;

                    Atom *atom2 = &*it2;
                    int index2 = get_index(atom2);

               if (j_atom <= i_atom) continue;
               if (chain_distance<ChainFB>(atom1,atom2) < 3) continue;
               if (atom2->mass == definitions::atom_h_weight) continue;

               if ((atom1->position - atom2->position).norm() > 9) continue;
                    const double energy_sum_temp = calculate_contribution(atom1, atom2,index1,index2);
                    energy_sum += energy_sum_temp;

                     //printf("ASC: EEF1-SB %4d %4d  Etot = %15.10f  Etemp = %15.10f  r_ij = %14.10f\n",
                     //        i_atom, j_atom, energy_sum, energy_sum_temp, (atom1->position - atom2->position).norm());
               }

          }

          printf("          EEF1-SB E = %15.6f kJ/mol\n", energy_sum * charmm36_constants::KCAL_TO_KJ);
          printf("          EEF1-SB E = %15.6f kcal/mol\n", energy_sum);
          return energy_sum;
     }

     };

}

#endif
