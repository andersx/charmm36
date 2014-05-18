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

#include <boost/tokenizer.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "charmm22_parser.h"

namespace phaistos {

const double exp_eef1 [350] = {
          0.99998,0.99978,0.99938,0.99878,0.99798,0.99698,
          0.99578,0.99439,0.99280,0.99102,0.98904,0.98686,0.98450,
          0.98194,0.97919,0.97626,0.97314,0.96984,0.96635,0.96269,
          0.95885,0.95483,0.95064,0.94627,0.94174,0.93704,0.93218,
          0.92716,0.92199,0.91665,0.91117,0.90554,0.89976,0.89384,
          0.88779,0.88159,0.87527,0.86882,0.86224,0.85554,0.84872,
          0.84179,0.83475,0.82760,0.82035,0.81300,0.80555,0.79802,
          0.79039,0.78268,0.77490,0.76703,0.75910,0.75109,0.74303,
          0.73490,0.72671,0.71847,0.71019,0.70186,0.69349,0.68508,
          0.67663,0.66816,0.65966,0.65114,0.64261,0.63405,0.62549,
          0.61691,0.60834,0.59976,0.59119,0.58262,0.57406,0.56551,
          0.55698,0.54847,0.53998,0.53151,0.52308,0.51467,0.50630,
          0.49797,0.48967,0.48142,0.47321,0.46504,0.45693,0.44887,
          0.44086,0.43291,0.42502,0.41719,0.40942,0.40171,0.39407,
          0.38650,0.37900,0.37157,0.36421,0.35693,0.34972,0.34259,
          0.33554,0.32856,0.32167,0.31486,0.30813,0.30149,0.29493,
          0.28845,0.28206,0.27576,0.26954,0.26341,0.25737,0.25142,
          0.24556,0.23978,0.23410,0.22850,0.22299,0.21757,0.21224,
          0.20700,0.20185,0.19679,0.19181,0.18693,0.18213,0.17742,
          0.17280,0.16826,0.16381,0.15945,0.15517,0.15098,0.14687,
          0.14284,0.13890,0.13503,0.13125,0.12755,0.12393,0.12039,
          0.11692,0.11354,0.11023,0.10699,0.10383,0.10074,0.09772,
          0.09478,0.09190,0.08910,0.08636,0.08369,0.08109,0.07855,
          0.07608,0.07367,0.07132,0.06903,0.06680,0.06463,0.06252,
          0.06047,0.05847,0.05653,0.05464,0.05280,0.05102,0.04928,
          0.04760,0.04596,0.04437,0.04283,0.04133,0.03987,0.03846,
          0.03710,0.03577,0.03449,0.03324,0.03203,0.03086,0.02973,
          0.02863,0.02757,0.02654,0.02555,0.02458,0.02365,0.02275,
          0.02188,0.02104,0.02023,0.01944,0.01869,0.01795,0.01725,
          0.01656,0.01590,0.01527,0.01465,0.01406,0.01349,0.01294,
          0.01241,0.01190,0.01141,0.01094,0.01048,0.01004,0.00962,
          0.00921,0.00882,0.00844,0.00808,0.00773,0.00740,0.00708,
          0.00677,0.00647,0.00619,0.00592,0.00565,0.00540,0.00516,
          0.00493,0.00470,0.00449,0.00429,0.00409,0.00390,0.00372,
          0.00355,0.00339,0.00323,0.00308,0.00293,0.00279,0.00266,
          0.00253,0.00241,0.00230,0.00219,0.00208,0.00198,0.00188,
          0.00179,0.00170,0.00162,0.00154,0.00146,0.00139,0.00132,
          0.00125,0.00119,0.00113,0.00107,0.00102,0.00097,0.00092,
          0.00087,0.00082,0.00078,0.00074,0.00070,0.00066,0.00063,
          0.00060,0.00056,0.00053,0.00051,0.00048,0.00045,0.00043,
          0.00040,0.00038,0.00036,0.00034,0.00032,0.00031,0.00029,
          0.00027,0.00026,0.00024,0.00023,0.00022,0.00020,0.00019,
          0.00018,0.00017,0.00016,0.00015,0.00014,0.00014,0.00013,
          0.00012,0.00011,0.00011,0.00010,0.00009,0.00009,0.00008,
          0.00008,0.00007,0.00007,0.00007,0.00006,0.00006,0.00005,
          0.00005,0.00005,0.00004,0.00004,0.00004,0.00003,0.00003,
          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
          0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
          0.0};

//! Gromacs van der Waals interaction term
class TermGromacsVdw: public EnergyTermCommon<TermGromacsVdw, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermGromacsVdw, ChainFB> EnergyTermCommon;

     //! Number of interactions calculated
     int counter;

     //! Lookup tables containing parameters
     std::vector<double> dGref;
     std::vector< std::vector<double> > factors;
     std::vector<double> vdw_radii;
     std::vector<double> lambda;

     std::map<std::string, unsigned int> eef1_atom_type_index_map;

     double dGref_total;
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


              initialize();
              std::string non_bonded_filename = "/home/andersx/phaistos_dev/modules/gromacs/src/energy/charmm22_cmap/charmm22_vdw.itp";
              non_bonded_parameters = read_nonbonded_parameters(non_bonded_filename);

              std::string non_bonded_14_filename = "/home/andersx/phaistos_dev/modules/gromacs/src/energy/charmm22_cmap/charmm22_vdw14.itp";
              non_bonded_14_parameters = read_nonbonded_14_parameters(non_bonded_14_filename);

              // std::cout << non_bonded_parameters[1].atom_type << std::endl;
              // std::cout << non_bonded_14_parameters[1].atom_type1 << std::endl;

              non_bonded_pairs = generate_non_bonded_pairs(this->chain,
                                                           non_bonded_parameters,
                                                           non_bonded_14_parameters,
                                                           dGref,
                                                           factors,
                                                           vdw_radii,
                                                           lambda,
                                                           eef1_atom_type_index_map);

            for (AtomIterator<ChainFB, definitions::ALL> it(*this->chain); !it.end(); ++it) {

                Atom *atom = &*it;
                std::string atom_type = get_charmm36_atom_type(atom);

                unsigned int index = eef1_atom_type_index_map[atom_type];

                dGref_total += dGref[index];
            }
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

     void initialize() {

          // Read parameter file
          std::vector<std::string> atoms;
          std::vector< std::vector<double> > params;
          std::string elem_str;

          //! Useful constant
          const double two_pi_3_2 = 2.0*M_PI*sqrt(M_PI);
          const double phys_t = 298.15;

          std::string settings_solvation_filename = "/home/andersx/phaistos_dev/modules/gromacs/src/energy/solvpar_17.inp";
          std::ifstream f_h((settings_solvation_filename).c_str());
          if (!f_h) {
               printf("# ERROR: EEF1 input: unable to open file '%s'\n", (settings_solvation_filename).c_str() );
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

          // double vdw_energy = 0.0;
          // double coul_energy = 0.0;


          // double vdw14_energy = 0.0;
          // double coul14_energy = 0.0;


          double eef1_sb_energy = dGref_total;
          //  int counter = 0;

          // std::cout << non_bonded_pairs.size() << " terms total" << std::endl;
          #pragma omp parallel for reduction(+:energy_sum,eef1_sb_energy) schedule(static)
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

               energy_sum += vdw_energy_temp + coul_energy_temp;

               if ((pair.do_eef1) && (r_sq < 81.0)) {
               // if (pair.do_eef1) std::cout << "DO EEF1" << std::endl;

                    const double r_ij = std::sqrt(r_sq);

                    double R_min_i = pair.R_vdw_1;
                    double R_min_j = pair.R_vdw_2;

                    double lambda_i = pair.lambda1;
                    double lambda_j = pair.lambda2;

                    const double arg_ij = std::fabs((r_ij - R_min_i)/lambda_i);
                    const double arg_ji = std::fabs((r_ij - R_min_j)/lambda_j);

                    int bin_ij = int(arg_ij*100);
                    int bin_ji = int(arg_ji*100);

                    double exp_ij = 0.0;
                    double exp_ji = 0.0;

                    if (bin_ij < 350) exp_ij = exp_eef1[bin_ij];
                    if (bin_ji < 350) exp_ji = exp_eef1[bin_ji];

                    double cont_ij = -pair.fac_12*exp_ij/r_sq;
                    double cont_ji = -pair.fac_21*exp_ji/r_sq;

                    eef1_sb_energy += cont_ij + cont_ji;
                   //  std::cout << pair.atom1 << "\n";
                   //  std::cout << pair.atom2 << "\n";
                   // std::cout << "EEF1:" << pair.lambda1 << "\n";
                   // std::cout << "EEF1:" << pair.lambda2 << "\n";
                   // std::cout << "EEF1:" << pair.fac_12 << "\n";
                   // std::cout << "EEF1:" << pair.fac_21 << "\n";
                   // std::cout << "EEF1:" << pair.R_vdw_1 << "\n";
                   // std::cout << "EEF1:" << pair.R_vdw_2 << "\n";
                   // std::cout << "EEF1:" << cont_ij << "\n";
                   // std::cout << "EEF1:" << cont_ji << "\n";

               }
               // if (!pair.is_14_interaction) {
               //     coul_energy += coul_energy_temp;
               //     vdw_energy += vdw_energy_temp;
               //  } else {
               //     coul14_energy += coul_energy_temp;
               //     vdw14_energy += vdw_energy_temp;

               //     // ::cout << counter++ << std::endl;

               //  // printf("ASC: LxCOL: %5d %5d %5d  r = %8.4f  XYZ = %8.4f %8.4f %8.4f   XYZ2 = %8.4f %8.4f %8.4f   q1 = %6.3f  q2 = %6.3f  qq = %7.3f  ecoul = %7.3f\n", pair.i1, pair.i2, pair.i2 , sqrt(r_sq),

               //     // printf("ASC: LxCOL: %5d  r = %8.4f  XYZ = %8.4f %8.4f %8.4f   XYZ2 = %8.4f %8.4f %8.4f   q1 = %6.3f  q2 = %6.3f  ecol = %7.3f  c6 = %24.20f  c12 = %24.20f  evdw = %14.10f\n",
               //     //      counter,
               //     //      sqrt(r_sq),
               //     //    (pair.atom1)->position[0],  (pair.atom1)->position[1],  (pair.atom1)->position[2],
               //     //    (pair.atom2)->position[0],  (pair.atom2)->position[1],  (pair.atom2)->position[2],
               //     //   pair.q1, pair.q2, coul_energy_temp, pair.c6, pair.c12, vdw_energy_temp);

               //     // counter++;
               //    //std::cout << "ASC: " << pair.i1
               //    // << "   " << pair.i2
               //    // // << "   " << pair.atom1
               //    // // << "   " << pair.atom2
               //    // << "  q1 = " << pair.q1
               //    // << "  q2 = " << pair.q2
               //    // << "  qq = " << pair.qq
               //    // // << "  s1 = " << pair.sigma1
               //    // // << "  s2 = " << pair.sigma2
               //    // // << "  e1 = " << pair.epsilon1
               //    // // << "  e2 = " << pair.epsilon2
               //    // << "   r = " << 10.0 / sqrt(inv_r_sq)
               //    // // << std::setprecision(9) << "   c6 = "   << pair.c6
               //    // // << "   c12 = "  << pair.c12
               //    // // << "   evdw = "  << vdw_energy_temp
               //    // << "   ecoul = "  << coul_energy_temp
               //    // // << "   TOTAL = "  << vdw_energy
               //    // << std::endl;

               //  }
          }



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

          // std::cout << "          coulomb E = " << coul_energy <<  " kcal/mol" << std::endl;
          // std::cout << "              vdW E = " << vdw_energy <<  " kcal/mol" << std::endl;
          // printf("       coulomb-14 E = %12.4f kJ/mol\n", coul14_energy);
          // printf("           vdW-14 E = %12.4f kJ/mol\n", vdw14_energy);
          // printf("       coulomb-SR E = %12.4f kJ/mol\n", coul_energy);
          // printf("           vdW-SR E = %12.4f kJ/mol\n", vdw_energy);
          // printf("          EEF1-SB E = %12.4f kcal/mol\n", eef1_sb_energy);
          // GROMACS energies are in kJ/mol, and need to return in kcal/mol
          //const double total_energy_in_kcal_mol = (coul_energy + vdw_energy) / 4.184;
          //return total_energy_in_kcal_mol;

          return energy_sum / 4.184 + eef1_sb_energy;
     }
};

}

#endif
