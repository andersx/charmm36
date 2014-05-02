// term_eef1.h ---  eef1-like solvation energy term
// Copyright (C) 2009-2011 Sandro Bottaro
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

#ifndef TERM_EEF1_H
#define TERM_EEF1_H

#include <boost/type_traits/is_base_of.hpp>
#include <boost/tokenizer.hpp>
#include "energy/energy_term.h"


int exponential [350] = { \
          0.99998,0.99978,0.99938,0.99878,0.99798,0.99698, \
          0.99578,0.99439,0.99280,0.99102,0.98904,0.98686,0.98450, \
          0.98194,0.97919,0.97626,0.97314,0.96984,0.96635,0.96269, \
          0.95885,0.95483,0.95064,0.94627,0.94174,0.93704,0.93218, \
          0.92716,0.92199,0.91665,0.91117,0.90554,0.89976,0.89384, \
          0.88779,0.88159,0.87527,0.86882,0.86224,0.85554,0.84872, \
          0.84179,0.83475,0.82760,0.82035,0.81300,0.80555,0.79802, \
          0.79039,0.78268,0.77490,0.76703,0.75910,0.75109,0.74303, \
          0.73490,0.72671,0.71847,0.71019,0.70186,0.69349,0.68508, \
          0.67663,0.66816,0.65966,0.65114,0.64261,0.63405,0.62549, \
          0.61691,0.60834,0.59976,0.59119,0.58262,0.57406,0.56551, \
          0.55698,0.54847,0.53998,0.53151,0.52308,0.51467,0.50630, \
          0.49797,0.48967,0.48142,0.47321,0.46504,0.45693,0.44887, \
          0.44086,0.43291,0.42502,0.41719,0.40942,0.40171,0.39407, \
          0.38650,0.37900,0.37157,0.36421,0.35693,0.34972,0.34259, \
          0.33554,0.32856,0.32167,0.31486,0.30813,0.30149,0.29493, \
          0.28845,0.28206,0.27576,0.26954,0.26341,0.25737,0.25142, \
          0.24556,0.23978,0.23410,0.22850,0.22299,0.21757,0.21224, \
          0.20700,0.20185,0.19679,0.19181,0.18693,0.18213,0.17742, \
          0.17280,0.16826,0.16381,0.15945,0.15517,0.15098,0.14687, \
          0.14284,0.13890,0.13503,0.13125,0.12755,0.12393,0.12039, \
          0.11692,0.11354,0.11023,0.10699,0.10383,0.10074,0.09772, \
          0.09478,0.09190,0.08910,0.08636,0.08369,0.08109,0.07855, \
          0.07608,0.07367,0.07132,0.06903,0.06680,0.06463,0.06252, \
          0.06047,0.05847,0.05653,0.05464,0.05280,0.05102,0.04928, \
          0.04760,0.04596,0.04437,0.04283,0.04133,0.03987,0.03846, \
          0.03710,0.03577,0.03449,0.03324,0.03203,0.03086,0.02973, \
          0.02863,0.02757,0.02654,0.02555,0.02458,0.02365,0.02275, \
          0.02188,0.02104,0.02023,0.01944,0.01869,0.01795,0.01725, \
          0.01656,0.01590,0.01527,0.01465,0.01406,0.01349,0.01294, \
          0.01241,0.01190,0.01141,0.01094,0.01048,0.01004,0.00962, \
          0.00921,0.00882,0.00844,0.00808,0.00773,0.00740,0.00708, \
          0.00677,0.00647,0.00619,0.00592,0.00565,0.00540,0.00516, \
          0.00493,0.00470,0.00449,0.00429,0.00409,0.00390,0.00372, \
          0.00355,0.00339,0.00323,0.00308,0.00293,0.00279,0.00266, \
          0.00253,0.00241,0.00230,0.00219,0.00208,0.00198,0.00188, \
          0.00179,0.00170,0.00162,0.00154,0.00146,0.00139,0.00132, \
          0.00125,0.00119,0.00113,0.00107,0.00102,0.00097,0.00092, \
          0.00087,0.00082,0.00078,0.00074,0.00070,0.00066,0.00063, \
          0.00060,0.00056,0.00053,0.00051,0.00048,0.00045,0.00043, \
          0.00040,0.00038,0.00036,0.00034,0.00032,0.00031,0.00029, \
          0.00027,0.00026,0.00024,0.00023,0.00022,0.00020,0.00019, \
          0.00018,0.00017,0.00016,0.00015,0.00014,0.00014,0.00013, \
          0.00012,0.00011,0.00011,0.00010,0.00009,0.00009,0.00008, \
          0.00008,0.00007,0.00007,0.00007,0.00006,0.00006,0.00005, \
          0.00005,0.00005,0.00004,0.00004,0.00004,0.00003,0.00003, \
          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, \
          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, \
          0.0};

namespace phaistos {

//! partial eef1 interaction term
class TermEef1: public EnergyTermCommon<TermEef1, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermEef1, ChainFB> EnergyTermCommon;

     //! Number of interactions in the last evaluation
     int counter;

protected:

     //! Lookup tables containing parameters
     std::vector<double> dGref;
     std::vector< std::vector<double> > factors;
     std::vector<double> vdw_radii;
     std::vector<double> lambda;

public:


     //! Local settings class
     const class Settings: public EnergyTerm<ChainFB>::SettingsClassicEnergy {
     public:

          //! Path to file containing solvation parameters
          std::string solvation_filename;

          //! Temperature
          double temp;

          //! Constructor
          Settings(std::string solvation_filename="/home/sandro/software/phaistos_1/phaistos/modules/force_field/data/solvpar_17.inp",double temp=298.15)
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
     TermEef1(ChainFB *chain,
              const Settings &settings=Settings(),
              RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "eef1", settings, random_number_engine) {

          initialize();

     }


     //! Copy constructor
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermEef1(const TermEef1 &other,
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
               if(abs(dt)>0.0001){
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
               if(abs(params[i][1])>0.001){
                    dGfree_i = (dGref[i]/(params[i][1]))*(params[i][2]);
               }

               std::vector<double> factors_i;
               for (unsigned int j=0; j<atoms.size(); j++){

                    //Factor = dGfree_i*V_j/(2*pi*sqrt(pi)*lambda_i)
                    double factor = 0.0;

                    if(abs(params[i][5])>0.001){
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

          int d = chain_distance<ChainFB>(atom1,atom2); //6s @ 500 iterations
          if (d > 3) {
               return calc_eef1_energy(atom1,atom2,index1,index2);
          } else {
               return 0.0;
          }
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

          double arg_ij = abs((r_ij - R_min_i)/lambda_i);
          double arg_ji = abs((r_ij - R_min_j)/lambda_j);
//          std::cout<< "rij " << r_ij << "\n";
//          std::cout<< "rijsq " << r_ij_sq << "\n";
//          std::cout<< "RMini " << R_min_i << "\n";
//          std::cout<< "RMinj " << R_min_j << "\n";
//          std::cout<< "lambdi " << lambda_i << "\n";
//          std::cout<< "lambdj " << lambda_j << "\n";
//
          // In CHARMM the exponential is not calculated explicitly
          // A lookup table is used instead (copied here for testing purposes only)
          int bin_ij = int(arg_ij*100);
          int bin_ji = int(arg_ji*100);
          double exp_ij = exponential[bin_ij];
          double exp_ji = exponential[bin_ji];
//          std::cout<< "expij " << exp_ij << "\n";
//          std::cout<< "expji " << exp_ji << "\n";
//          std::cout<< "facij " << R_min_i << "\n";
//          std::cout<< "facji " << factors[index1][index2] << "\n";
//          std::cout<< "facji " << factors[index2][index1] << "\n";

          double cont_ij = -factors[index1][index2]*exp_ij/r_ij_sq;
          double cont_ji = -factors[index2][index1]*exp_ji/r_ij_sq;

          // Use standard exponentiation (no lookup table)
          //double cont_ij = -factors[index1][index2]*exp(-(arg_ij*arg_ij));
          //double cont_ji = -factors[index1][index2]*exp(-(arg_ji*arg_ji));

          return (cont_ij+cont_ji);
     }

     //! Return index in parameter table corresponding to the atomtype
     int get_index(Atom *atom){
          // Return 5 now (5 is nice) - To be corrected!
          return 20;
     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return vdw potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energy_sum=0.0;
          counter=0;

          // Iterate all the atom pairs on the chain
          for (AtomIterator<ChainFB, definitions::ALL> it1(*this->chain); !it1.end(); ++it1) {

               // Add dGref contribution
               // Remember to skip hydrogens (and atoms with dGref == 0)

               Atom *atom1 = &*it1;
               int index1 = get_index(atom1);
               energy_sum += dGref[index1];

               for(AtomIterator<ChainFB, definitions::ALL> it2(it1+1); !it2.end(); ++it2){

                    Atom *atom2 = &*it2;
                    int index2 = get_index(atom2);
                    energy_sum += calculate_contribution(atom1, atom2,index1,index2);
               }

          }
          return energy_sum;
     }

     };

}

#endif
