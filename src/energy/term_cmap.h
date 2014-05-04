// term_cmap.h --- CMAP torsion angle energy term
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

#ifndef TERM_GROMACS_CMAP_H
#define TERM_GROMACS_CMAP_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"


namespace gromacs_cmap {

const int cmap_coeff_matrix[] = {
    1, 0, -3,  2, 0, 0,  0,  0, -3,  0,  9, -6,  2,  0, -6,  4,
    0, 0,  0,  0, 0, 0,  0,  0,  3,  0, -9,  6, -2,  0,  6, -4,
    0, 0,  0,  0, 0, 0,  0,  0,  0,  0,  9, -6,  0,  0, -6,  4,
    0, 0,  3, -2, 0, 0,  0,  0,  0,  0, -9,  6,  0,  0,  6, -4,
    0, 0,  0,  0, 1, 0, -3,  2, -2,  0,  6, -4,  1,  0, -3,  2,
    0, 0,  0,  0, 0, 0,  0,  0, -1,  0,  3, -2,  1,  0, -3,  2,
    0, 0,  0,  0, 0, 0,  0,  0,  0,  0, -3,  2,  0,  0,  3, -2,
    0, 0,  0,  0, 0, 0,  3, -2,  0,  0, -6,  4,  0,  0,  3, -2,
    0, 1, -2,  1, 0, 0,  0,  0,  0, -3,  6, -3,  0,  2, -4,  2,
    0, 0,  0,  0, 0, 0,  0,  0,  0,  3, -6,  3,  0, -2,  4, -2,
    0, 0,  0,  0, 0, 0,  0,  0,  0,  0, -3,  3,  0,  0,  2, -2,
    0, 0, -1,  1, 0, 0,  0,  0,  0,  0,  3, -3,  0,  0, -2,  2,
    0, 0,  0,  0, 0, 1, -2,  1,  0, -2,  4, -2,  0,  1, -2,  1,
    0, 0,  0,  0, 0, 0,  0,  0,  0, -1,  2, -1,  0,  1, -2,  1,
    0, 0,  0,  0, 0, 0,  0,  0,  0,  0,  1, -1,  0,  0, -1,  1,
    0, 0,  0,  0, 0, 0, -1,  1,  0,  0,  2, -2,  0,  0, -1,  1
};

const int loop_index[4][4] = {
        {0, 4, 8, 12},
        {1, 5, 9, 13},
        {2, 6, 10, 14},
        {3, 7, 11, 15}
    };

class cmap_dihedral_type() {

}


void parse_cmap_parameters()



} // end namespace gromacs_cmap

namespace phaistos {

//! Torsion energy term
class TermGromacsCmap: public EnergyTermCommon<TermGromacsTorsion, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermGromacsCmap, ChainFB> EnergyTermCommon;

     //! Number of interactions calculated
     int counter;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermGromacsCmap(ChainFB *chain,
                        const Settings &settings=Settings(),
                        RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "gromacs-cmap", settings, random_number_engine) {

     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermGromacsCmap(const TermGromacsCmap &other,
                        RandomNumberEngine *random_number_engine,
                        int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            counter(other.counter) {

     }


    int cmap_setup_grid_index(int ip, int grid_spacing, int *ipm1, int *ipp1, int *ipp2) {

        int im1, ip1, ip2;

        if (ip < 0) {
            ip = ip + grid_spacing - 1;
        } else if (ip > grid_spacing) {
            ip = ip - grid_spacing - 1;
        }

        im1 = ip - 1;
        ip1 = ip + 1;
        ip2 = ip + 2;

        if (ip == 0){
            im1 = grid_spacing - 1;
        } else if (ip == grid_spacing-2) {
            ip2 = 0;
        } else if (ip == grid_spacing-1) {
            ip1 = 0;
            ip2 = 1;
        }

        *ipm1 = im1;
        *ipp1 = ip1;
        *ipp2 = ip2;

        return ip;

    }



    // CMAP table lookup and interpolation algorightm from GROMACS
    // Needs input:
    // Phi, Psi, dihedral type, parsed cmap.itp, gridspace
    double cmap_dihedral_energy() {

        double e  = 0;
        // int         i, j, k, n, idx;
        // int         ai, aj, ak, al, am;
        // int         a1i, a1j, a1k, a1l, a2i, a2j, a2k, a2l;
        // int         type, cmapA;
        // int         t11, t21, t31, t12, t22, t32;
        // int         iphi1, ip1m1, ip1p1, ip1p2;
        // int         iphi2, ip2m1, ip2p1, ip2p2;
        // int         l1, l2, l3, l4;
        // int         pos1, pos2, pos3, pos4, tmp;

        // real        ty[4], ty1[4], ty2[4], ty12[4], tc[16], tx[16];
        // real        phi1, psi1, cos_phi1, sin_phi1, sign1, xphi1;
        // real        phi2, psi2, cos_phi2, sin_phi2, sign2, xphi2;
        // real        dx, xx, tt, tu, e, df1, df2, ddf1, ddf2, ddf12, vtot;
        // real        ra21, rb21, rg21, rg1, rgr1, ra2r1, rb2r1, rabr1;
        // real        ra22, rb22, rg22, rg2, rgr2, ra2r2, rb2r2, rabr2;
        // real        fg1, hg1, fga1, hgb1, gaa1, gbb1;
        // real        fg2, hg2, fga2, hgb2, gaa2, gbb2;
        // real        fac;

        // rvec        r1_ij, r1_kj, r1_kl, m1, n1;
        // rvec        r2_ij, r2_kj, r2_kl, m2, n2;
        // rvec        f1_i, f1_j, f1_k, f1_l;
        // rvec        f2_i, f2_j, f2_k, f2_l;
        // rvec        a1, b1, a2, b2;
        // rvec        f1, g1, h1, f2, g2, h2;
        // rvec        dtf1, dtg1, dth1, dtf2, dtg2, dth2;
        // ivec        jt1, dt1_ij, dt1_kj, dt1_lj;
        // ivec        jt2, dt2_ij, dt2_kj, dt2_lj;

        // The cmap lookup table
        // const real *cmapd;

        // xphi1, xphi2 // = Phi, psi in radians


        // // Number of grid pointdds
        // dx = 2*M_PI / cmap_grid->grid_spacing;

        // // Where on the grid are we
        // iphi1 = (int)(xphi1/dx);
        // iphi2 = (int)(xphi2/dx);

        // iphi1 = cmap_setup_grid_index(iphi1, cmap_grid->grid_spacing, &ip1m1, &ip1p1, &ip1p2);
        // iphi2 = cmap_setup_grid_index(iphi2, cmap_grid->grid_spacing, &ip2m1, &ip2p1, &ip2p2);

        // pos1    = iphi1*cmap_grid->grid_spacing+iphi2;
        // pos2    = ip1p1*cmap_grid->grid_spacing+iphi2;
        // pos3    = ip1p1*cmap_grid->grid_spacing+ip2p1;
        // pos4    = iphi1*cmap_grid->grid_spacing+ip2p1;

        // ty[0]   = cmapd[pos1*4];
        // ty[1]   = cmapd[pos2*4];
        // ty[2]   = cmapd[pos3*4];
        // ty[3]   = cmapd[pos4*4];

        // ty1[0]   = cmapd[pos1*4+1];
        // ty1[1]   = cmapd[pos2*4+1];
        // ty1[2]   = cmapd[pos3*4+1];
        // ty1[3]   = cmapd[pos4*4+1];

        // ty2[0]   = cmapd[pos1*4+2];
        // ty2[1]   = cmapd[pos2*4+2];
        // ty2[2]   = cmapd[pos3*4+2];
        // ty2[3]   = cmapd[pos4*4+2];

        // ty12[0]   = cmapd[pos1*4+3];
        // ty12[1]   = cmapd[pos2*4+3];
        // ty12[2]   = cmapd[pos3*4+3];
        // ty12[3]   = cmapd[pos4*4+3];

        // dx    = 360.0 / cmap_grid->grid_spacing;
        // xphi1 = xphi1 * RAD2DEG;
        // xphi2 = xphi2 * RAD2DEG;

        // for (i = 0; i < 4; i++) {
        //     tx[i]    = ty[i];
        //     tx[i+4]  = ty1[i]*dx;
        //     tx[i+8]  = ty2[i]*dx;
        //     tx[i+12] = ty12[i]*dx*dx;
        // }

        // idx = 0;
        // for (i = 0; i < 4; i++) {
        //     for (j = 0; j < 4; j++) {
        //         xx = 0;
        //         for (k = 0; k < 16; k++) {
        //             xx = xx + cmap_coeff_matrix[k*16+idx]*tx[k];
        //         }

        //         idx++;
        //         tc[i*4+j] = xx;
        //     }
        // }

        // tt = (xphi1-iphi1*dx)/dx;
        // tu = (xphi2-iphi2*dx)/dx;


        // for (unsigned int i = 3; i >= 0; i--) {

        //     l1 = loop_index[i][3];
        //     l2 = loop_index[i][2];
        //     l3 = loop_index[i][1];
        //     e  = tt * e  + ((tc[i*4+3]*tu+tc[i*4+2])*tu + tc[i*4+1])*tu+tc[i*4];
        // }

        return e;

    }


     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energy = 0.0;

          for (AtomIterator<ChainFB, definitions::ALL> it1(*this->chain);
               !it1.end(); ++it1) {

          }



          return energy;
     }
};

}

#endif
