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

#ifndef TERM_CHARMM36_CMAP_H
#define TERM_CHARMM36_CMAP_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "protein/iterators/pair_iterator_chaintree.h"
#include "protein/chain_fb.h"
#include "protein/definitions.h"

#include "term_cmap_tables.h"

namespace phaistos {

//! CMAP energy term
class TermCharmm36Cmap: public EnergyTermCommon<TermCharmm36Cmap, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36Cmap, ChainFB> EnergyTermCommon;

public:

     //double cmapdata[6][2304];
     std::vector<std::vector<double> > cmapdata;

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36Cmap(ChainFB *chain,
                        const Settings &settings=Settings(),
                        RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-cmap", settings, random_number_engine) {
          setup_cmap();
          }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36Cmap(const TermCharmm36Cmap &other,
                        RandomNumberEngine *random_number_engine,
                        int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            cmapdata(other.cmapdata) {}


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


void
spline1d( double        dx,
		 double *      y,
		 int           n,
		 double *      u,
		 double *      y2 )
{
    int i;
    double p,q;
	
    y2[0] = 0.0;
    u[0]  = 0.0;
	
    for(i=1;i<n-1;i++)
    {
		p = 0.5*y2[i-1]+2.0;
        y2[i] = -0.5/p;
        q = (y[i+1]-2.0*y[i]+y[i-1])/dx;
		u[i] = (3.0*q/dx-0.5*u[i-1])/p;
    }
	
    y2[n-1] = 0.0;
	
    for(i=n-2;i>=0;i--)
    {
        y2[i] = y2[i]*y2[i+1]+u[i];
    }
}


void
interpolate1d( double     xmin,
			  double     dx,
			  double *   ya,
			  double *   y2a,
			  double     x,
			  double *   y,
			  double *   y1)
{
    int ix;
    double a,b;
	
    ix = (x-xmin)/dx;
	
    a = (xmin+(ix+1)*dx-x)/dx;
    b = (x-xmin-ix*dx)/dx;
	
    *y  = a*ya[ix]+b*ya[ix+1]+((a*a*a-a)*y2a[ix]+(b*b*b-b)*y2a[ix+1])*(dx*dx)/6.0;
    *y1 = (ya[ix+1]-ya[ix])/dx-(3.0*a*a-1.0)/6.0*dx*y2a[ix]+(3.0*b*b-1.0)/6.0*dx*y2a[ix+1];
}


void setup_cmap () {

    int nc = 6;
    int grid_spacing = 24;

    const double * grid = charmm36_cmap::cmapd;

    this->cmapdata.resize(nc);
    for (int i = 0; i < nc; ++i)
            this->cmapdata[i].resize(grid_spacing*grid_spacing*4);
    printf("CMAP SETUP %5d  %f\n", nc, grid[1]);

 	// double *tmp_u,*tmp_u2,*tmp_yy,*tmp_y1,*tmp_t2,*tmp_grid;
	
    int    i,j,k,ii,jj,kk,idx;
	int    offset;
    double dx,xmin,v,v1,v2,v12;
    double phi,psi;
	
	// snew(tmp_u,2*grid_spacing);
	// snew(tmp_u2,2*grid_spacing);
	// snew(tmp_yy,2*grid_spacing);
	// snew(tmp_y1,2*grid_spacing);
	// snew(tmp_t2,2*grid_spacing*2*grid_spacing);
	// snew(tmp_grid,2*grid_spacing*2*grid_spacing);

    double *tmp_u       = new double[2*grid_spacing];
	double *tmp_u2      = new double[2*grid_spacing];
	double *tmp_yy      = new double[2*grid_spacing];
	double *tmp_y1      = new double[2*grid_spacing];
	double *tmp_t2      = new double[2*grid_spacing*2*grid_spacing];
	double *tmp_grid    = new double[2*grid_spacing*2*grid_spacing];

    dx = 360.0/grid_spacing;
    xmin = -180.0-dx*grid_spacing/2;
	
	for(kk=0;kk<nc;kk++)
	{
		/* Compute an offset depending on which cmap we are using                                 
		 * Offset will be the map number multiplied with the grid_spacing * grid_spacing * 2      
		 */
		offset = kk * grid_spacing * grid_spacing * 2;
		
		for(i=0;i<2*grid_spacing;i++)
		{
			ii=(i+grid_spacing-grid_spacing/2)%grid_spacing;
			
			for(j=0;j<2*grid_spacing;j++)
			{
				jj=(j+grid_spacing-grid_spacing/2)%grid_spacing;
				tmp_grid[i*grid_spacing*2+j] = grid[offset+ii*grid_spacing+jj];
			}
		}
		
		for(i=0;i<2*grid_spacing;i++)
		{
			spline1d(dx,&(tmp_grid[2*grid_spacing*i]),2*grid_spacing,tmp_u,&(tmp_t2[2*grid_spacing*i]));
		}
		
		for(i=grid_spacing/2;i<grid_spacing+grid_spacing/2;i++)
		{
			ii = i-grid_spacing/2;
			phi = ii*dx-180.0;
			
			for(j=grid_spacing/2;j<grid_spacing+grid_spacing/2;j++)
			{
				jj = j-grid_spacing/2;
				psi = jj*dx-180.0;
				
				for(k=0;k<2*grid_spacing;k++)
				{
					interpolate1d(xmin,dx,&(tmp_grid[2*grid_spacing*k]),
								  &(tmp_t2[2*grid_spacing*k]),psi,&tmp_yy[k],&tmp_y1[k]);
				}
				
				spline1d(dx,tmp_yy,2*grid_spacing,tmp_u,tmp_u2);
				interpolate1d(xmin,dx,tmp_yy,tmp_u2,phi,&v,&v1);
				spline1d(dx,tmp_y1,2*grid_spacing,tmp_u,tmp_u2);
				interpolate1d(xmin,dx,tmp_y1,tmp_u2,phi,&v2,&v12);
				
				idx = ii*grid_spacing+jj;
				this->cmapdata[kk][idx*4] = grid[offset+ii*grid_spacing+jj];
				this->cmapdata[kk][idx*4+1] = v1;
				this->cmapdata[kk][idx*4+2] = v2;
				this->cmapdata[kk][idx*4+3] = v12;

                printf("cmap_grid->cmapdata[%d].cma[%d] = %24.18f   grid[offset+ii*grid_spacing+jj] = %d\n", kk, idx*4+1, grid[offset+ii*grid_spacing+jj], offset+ii*grid_spacing+jj);
                printf("cmap_grid->cmapdata[%d].cma[%d] = %24.18f\n", kk, idx*4+1, v1);
                printf("cmap_grid->cmapdata[%d].cma[%d] = %24.18f\n", kk, idx*4+2, v2);
                printf("cmap_grid->cmapdata[%d].cma[%d] = %24.18f\n", kk, idx*4+3, v12);

                printf("CMAP VAR:   N = %d   POS = %d  val = %24.18f\n", kk, ii*grid_spacing+jj, grid[offset+ii*grid_spacing+jj]);
			}
		}
	}
}				



    // CMAP table lookup and interpolation algorightm from CHARMM36
    // Needs input:
    // Phi, Psi, dihedral type, parsed cmap.itp, gridspace
    double cmap_dihedral_energy(const double phi, const double psi) {

         double e  = 0;
         // int         i, j, k, n, idx;
         int         i, j, k, idx;
         //int         ai, aj, ak, al, am;
         //int         a1i, a1j, a1k, a1l, a2i, a2j, a2k, a2l;
         //int         type, cmapA;
         //int         t11, t21, t31, t12, t22, t32;
         int         iphi1, ip1m1, ip1p1, ip1p2;
         int         iphi2, ip2m1, ip2p1, ip2p2;
         //int         l1, l2, l3, l4;
         int         pos1, pos2, pos3, pos4;// , tmp;

         double        ty[4], ty1[4], ty2[4], ty12[4], tc[16], tx[16];
         // double        phi1, psi1, cos_phi1, sin_phi1, sign1, xphi1;
         // double        phi2, psi2, cos_phi2, sin_phi2, sign2, xphi2;
         double        phi1,  xphi1;
         double        phi2,  xphi2;
         double        dx, xx, tt, tu; // df1, df2, ddf1, ddf2, ddf12, vtot;
         //double        ra21, rb21, rg21, rg1, rgr1, ra2r1, rb2r1, rabr1;
         //double        ra22, rb22, rg22, rg2, rgr2, ra2r2, rb2r2, rabr2;
         //double        fg1, hg1, fga1, hgb1, gaa1, gbb1;
         //double        fg2, hg2, fga2, hgb2, gaa2, gbb2;
         //double        fac;


         phi1 = phi;
         phi2 = psi;

         const int grid_spacing = 24;
         const double RAD2DEG = 180.0 / M_PI;

		 xphi1 = phi1 + M_PI;
		 xphi2 = phi2 + M_PI;

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

         // xphi1, xphi2 // = Phi, psi in radians


         // Number of grid pointdds
         dx = 2*M_PI / grid_spacing;

         printf("CMAP: GRID %5d\n", grid_spacing);

         // Where on the grid are we
         iphi1 = (int)(xphi1/dx);
         iphi2 = (int)(xphi2/dx);

         iphi1 = cmap_setup_grid_index(iphi1, grid_spacing, &ip1m1, &ip1p1, &ip1p2);
         iphi2 = cmap_setup_grid_index(iphi2, grid_spacing, &ip2m1, &ip2p1, &ip2p2);


         printf("\nASC: CMAP   iphi1 = %4d   iphi2 = %4d\n", iphi1, iphi2);
         pos1    = iphi1*grid_spacing+iphi2;
         pos2    = ip1p1*grid_spacing+iphi2;
         pos3    = ip1p1*grid_spacing+ip2p1;
         pos4    = iphi1*grid_spacing+ip2p1;


         // printf("CMAPPOS: %4d  %4d  %4d %4d %14.8f\n", pos1, pos2, pos3, pos4, charmm36_cmap::cmapd[2]);
         // ty[0]   = charmm36_cmap::cmapd[pos1*4];
         // ty[1]   = charmm36_cmap::cmapd[pos2*4];
         // ty[2]   = charmm36_cmap::cmapd[pos3*4];
         // ty[3]   = charmm36_cmap::cmapd[pos4*4];

         // printf("CMAPTY %14.8f  %14.8f  %14.8f  %14.8f\n", ty[0], ty[1], ty[2], ty[3]);
         // ty1[0]   = charmm36_cmap::cmapd[pos1*4+1];
         // ty1[1]   = charmm36_cmap::cmapd[pos2*4+1];
         // ty1[2]   = charmm36_cmap::cmapd[pos3*4+1];
         // ty1[3]   = charmm36_cmap::cmapd[pos4*4+1];

         // ty2[0]   = charmm36_cmap::cmapd[pos1*4+2];
         // ty2[1]   = charmm36_cmap::cmapd[pos2*4+2];
         // ty2[2]   = charmm36_cmap::cmapd[pos3*4+2];
         // ty2[3]   = charmm36_cmap::cmapd[pos4*4+2];

         // ty12[0]   = charmm36_cmap::cmapd[pos1*4+3];
         // ty12[1]   = charmm36_cmap::cmapd[pos2*4+3];
         // ty12[2]   = charmm36_cmap::cmapd[pos3*4+3];
         // ty12[3]   = charmm36_cmap::cmapd[pos4*4+3];
         ty[0]     = this->cmapdata[0][pos1*4];
         ty[1]     = this->cmapdata[0][pos2*4];
         ty[2]     = this->cmapdata[0][pos3*4];
         ty[3]     = this->cmapdata[0][pos4*4];

         ty1[0]    = this->cmapdata[0][pos1*4+1];
         ty1[1]    = this->cmapdata[0][pos2*4+1];
         ty1[2]    = this->cmapdata[0][pos3*4+1];
         ty1[3]    = this->cmapdata[0][pos4*4+1];

         ty2[0]    = this->cmapdata[0][pos1*4+2];
         ty2[1]    = this->cmapdata[0][pos2*4+2];
         ty2[2]    = this->cmapdata[0][pos3*4+2];
         ty2[3]    = this->cmapdata[0][pos4*4+2];

         ty12[0]   = this->cmapdata[0][pos1*4+3];
         ty12[1]   = this->cmapdata[0][pos2*4+3];
         ty12[2]   = this->cmapdata[0][pos3*4+3];
         ty12[3]   = this->cmapdata[0][pos4*4+3];

         printf("CMAPTY %14.8f  %14.8f  %14.8f  %14.8f\n", ty[0], ty[1], ty[2], ty[3]);

         dx    = 360.0 / grid_spacing;
         xphi1 = xphi1 * RAD2DEG;
         xphi2 = xphi2 * RAD2DEG;


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
                     xx = xx + charmm36_cmap::cmap_coeff_matrix[k*16+idx]*tx[k];
                 }

                 idx++;
                 tc[i*4+j] = xx;
             }
         }

         tt = (xphi1-iphi1*dx)/dx;
         tu = (xphi2-iphi2*dx)/dx;


         for (int i = 3; i >= 0; i--) {

             //l1 = charmm36_cmap::loop_index[i][3];
             //l2 = charmm36_cmap::loop_index[i][2];
             //l3 = charmm36_cmap::loop_index[i][1];
             e  = tt * e  + ((tc[i*4+3]*tu+tc[i*4+2])*tu + tc[i*4+1])*tu+tc[i*4];
            printf("CMAP: HERE1 %14.8f\n", e);
         }

         printf("\nASC: CMAP   phi1 = %14.8f   phi2 = %14.8f  e_cmap = %14.8f \n",
                phi/M_PI*180.0, psi/M_PI*180.0, e);

         return e;


    }


     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double energy = 0.0;

          for (ResidueIterator<ChainFB> res(*(chain)); !(res).end(); ++res) {

                energy += cmap_dihedral_energy(res->get_phi(), res->get_psi());

          }



          return energy;
     }
};

}

#endif
