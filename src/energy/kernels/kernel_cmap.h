// kernel_cmap.h --- CMAP torsion angle energy term
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

#ifndef KERNEL_CMAP_H
#define KERNEL_CMAP_H

namespace kernel {

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

//! CMAP table lookup and interpolation algorithm from Gromacs.
//! \param phi Phi angle of the residue
//! \param psi Psi angle of the residue
//! \param cmap_type_index Index {0 ... 5} for which CMAP table to use.
//! \param cmapdata Table which contains the CMAP correction tables in Gromacs' internal format.
double cmap_energy(const double phi,
                   const double psi,
                   const unsigned int cmap_type_index,
                   const std::vector<std::vector<double> > &cmapdata) {

     // The next three lines are added to in this Phaistos wrapper.
     double e  = 0.0;
     const int grid_spacing = 24;
     const double gromacs_RAD2DEG = 180.0 / M_PI;

     int         i, j, k, idx;
     int         iphi1, ip1m1, ip1p1, ip1p2;
     int         iphi2, ip2m1, ip2p1, ip2p2;
     int         pos1, pos2, pos3, pos4;

     double        ty[4], ty1[4], ty2[4], ty12[4], tc[16], tx[16];
     double        xphi1;
     double        xphi2;
     double        dx, xx, tt, tu;

	 xphi1 = phi + M_PI;
	 xphi2 = psi + M_PI;

     // Range mangling
	 if(xphi1<0)
	 {
	 	xphi1 = xphi1 + 2*M_PI;
	 }
	 else if(xphi1>=2*M_PI)
	 {
	 	xphi1 = xphi1 - 2*M_PI;
	 }

	 if(xphi2<0)
	 {
	 	xphi2 = xphi2 + 2*M_PI;
	 }
	 else if(xphi2>=2*M_PI)
	 {
		xphi2 = xphi2 - 2*M_PI;
	 }

     // Number of grid points
     dx = 2*M_PI / grid_spacing;

     // Where on the grid are we
     iphi1 = (int)(xphi1/dx);
     iphi2 = (int)(xphi2/dx);

     iphi1 = kernel::cmap_setup_grid_index(iphi1, grid_spacing, &ip1m1, &ip1p1, &ip1p2);
     iphi2 = kernel::cmap_setup_grid_index(iphi2, grid_spacing, &ip2m1, &ip2p1, &ip2p2);

     pos1    = iphi1*grid_spacing+iphi2;
     pos2    = ip1p1*grid_spacing+iphi2;
     pos3    = ip1p1*grid_spacing+ip2p1;
     pos4    = iphi1*grid_spacing+ip2p1;

     ty[0]     = cmapdata[cmap_type_index][pos1*4];
     ty[1]     = cmapdata[cmap_type_index][pos2*4];
     ty[2]     = cmapdata[cmap_type_index][pos3*4];
     ty[3]     = cmapdata[cmap_type_index][pos4*4];

     ty1[0]    = cmapdata[cmap_type_index][pos1*4+1];
     ty1[1]    = cmapdata[cmap_type_index][pos2*4+1];
     ty1[2]    = cmapdata[cmap_type_index][pos3*4+1];
     ty1[3]    = cmapdata[cmap_type_index][pos4*4+1];

     ty2[0]    = cmapdata[cmap_type_index][pos1*4+2];
     ty2[1]    = cmapdata[cmap_type_index][pos2*4+2];
     ty2[2]    = cmapdata[cmap_type_index][pos3*4+2];
     ty2[3]    = cmapdata[cmap_type_index][pos4*4+2];

     ty12[0]   = cmapdata[cmap_type_index][pos1*4+3];
     ty12[1]   = cmapdata[cmap_type_index][pos2*4+3];
     ty12[2]   = cmapdata[cmap_type_index][pos3*4+3];
     ty12[3]   = cmapdata[cmap_type_index][pos4*4+3];

     dx    = 360.0 / grid_spacing;
     xphi1 = xphi1 * gromacs_RAD2DEG;
     xphi2 = xphi2 * gromacs_RAD2DEG;

     for (i = 0; i < 4; i++) {
         tx[i]    = ty[i];
         tx[i+4]  = ty1[i]*dx;
         tx[i+8]  = ty2[i]*dx;
         tx[i+12] = ty12[i]*dx*dx;
     }

     idx = 0;
     for (i = 0; i < 4; i++) {
         for (j = 0; j < 4; j++) {
             xx = 0;
             for (k = 0; k < 16; k++) {
                 xx = xx + cmap_coeff_matrix[k*16+idx]*tx[k];
             }

             idx++;
             tc[i*4+j] = xx;
         }
     }

     tt = (xphi1-iphi1*dx)/dx;
     tu = (xphi2-iphi2*dx)/dx;


     for (int i = 3; i >= 0; i--) {
         e  = tt * e  + ((tc[i*4+3]*tu+tc[i*4+2])*tu + tc[i*4+1])*tu+tc[i*4];
     }

     // printf("ASC: CMAP   type = %2d  phi1 = %14.8f   phi2 = %14.8f  e_cmap = %14.8f \n",
     //       cmap_type_index, phi/M_PI*180.0, psi/M_PI*180.0, e);

     return e;
}


} // End namespace kernel


#endif
