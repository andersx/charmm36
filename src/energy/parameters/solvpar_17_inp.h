// solvpar_inp.h --- std::string'ed version of solv_par.inp
// Copyright (C) 2014 Sandro Bottaro, Anders S. Christensen
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

#ifndef TERM_CHARMM_SOLVPAR_INP_H
#define TERM_CHARMM_SOLVPAR_INP_H

namespace charmm_constants {

const static std::string solvpar_inp = "! Parameters for solvation free energy calculation in EEF1\n\
! Volume in A^3, energies in Kcal/mol, cp in cal/molK\n\
! Sigw (A) is the distance within which 84% of Gsolv arises\n\
!\n\
! EEF1SB: Sandro Bottaro, Kresten Lindorff-Larsen, Robert B. Best. \"Variational\n\
! optimization of an all-atom implicit solvent force field to match\n\
! explicit solvent experimental data\". J. Chem. Theor. Comput.\n\
! dx.doi.org/10.1021/ct400730n\n\
!\n\
! The van der Waal radii are taken from CHARMM force field for this implementation.\n\
!\n\
!\n\
!       Volume       Gref        Gfree       Href        CPref      Sigw        vdw Radius\n\
H        0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     0.22450\n\
HA       0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     1.32000\n\
HA1      0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     1.34000\n\
HA2      0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     1.34000\n\
HA3      0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     1.34000\n\
HB1      0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     1.32000\n\
HB2      0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     1.34000\n\
HC       0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     0.22450\n\
HE1      0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     1.25000\n\
HE2      0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     1.26000\n\
HP       0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     1.35820\n\
HR1      0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     0.90000\n\
HR2      0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     0.70000\n\
HR3      0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     1.46800\n\
HS       0.00000     0.00000     0.00000     0.00000     0.00000    3.50000     0.45000\n\
!\n\
C       14.72040     0.00000     0.00000     0.00000     0.00000    3.50000     2.00000\n\
CD      14.72040     0.00000     0.00000     0.00000     0.00000    3.50000     2.00000\n\
CT1     11.50710    -0.18700    -0.18700     0.87600     0.00000    3.50000     2.00000\n\
CT2     18.85010     0.37200     0.37200    -0.61000    18.60000    3.50000     2.01000\n\
CT2A    18.66580     0.37200     0.37200    -0.61000    18.60000    3.50000     2.01000\n\
CT3     27.94090     1.08900     1.08900    -1.77900    35.60000    3.50000     2.04000\n\
CPH1     5.27500     0.05700     0.08000    -0.97300     6.90000    3.50000     1.80000\n\
CPH2    11.79580     0.05700     0.08000    -0.97300     6.90000    3.50000     1.80000\n\
CPT      4.66880    -0.89000    -0.89000     2.22000     6.90000    3.50000     1.86000\n\
CY      10.50670    -0.89000    -0.89000     2.22000     6.90000    3.50000     1.99000\n\
CP1     25.45760    -0.18700    -0.18700     0.87600     0.00000    3.50000     2.27500\n\
CP2     19.88040     0.37200     0.37200    -0.61000    18.60000    3.50000     2.17500\n\
CP3     26.73090     0.37200     0.37200    -0.61000    18.60000    3.50000     2.17500\n\
CC      16.53920     0.00000     0.00000     0.00000     0.00000    3.50000     2.00000\n\
CAI     18.24850     0.05700     0.05700    -0.97300     6.90000    3.50000     1.99000\n\
CA      18.24850     0.05700     0.05700    -0.97300     6.90000    3.50000     1.99240\n\
!\n\
N        0.00000    -1.00000    -1.00000    -1.25000     8.80000    3.50000     1.85000\n\
NR1     15.27340    -5.95000    -5.95000    -9.05900    -8.80000    3.50000     1.85000\n\
NR2     15.11100    -3.82000    -3.82000    -4.65400    -8.80000    3.50000     1.85000\n\
NR3     15.07050    -5.95000    -5.95000    -9.05900    -8.80000    3.50000     1.85000\n\
NH1     10.19720    -5.95000    -5.95000    -9.05900    -8.80000    3.50000     1.85000\n\
NH2     18.18240    -5.95000    -5.95000    -9.05900    -8.80000    3.50000     1.85000\n\
NH3     18.81700   -20.00000   -20.00000   -25.00000   -18.00000    6.00000     1.85000\n\
NC2     18.21540   -10.00000   -10.00000   -12.00000    -7.00000    6.00000     1.85000\n\
NY      12.00100    -5.95000    -5.95000    -9.05900    -8.80000    3.50000     1.85000\n\
NP       4.99280   -20.00000   -20.00000   -25.00000   -18.00000    6.00000     1.85000\n\
!\n\
O       11.77220    -5.33000    -5.33000    -5.78700    -8.80000    3.50000     1.70000\n\
OB      11.69390    -5.33000    -5.33000    -5.78700    -8.80000    3.50000     1.70000\n\
OC      12.00340   -10.00000   -10.00000   -12.00000    -9.40000    6.00000     1.70000\n\
OH1     15.52750    -5.92000    -5.92000    -9.26400   -11.20000    3.50000     1.77000\n\
OS       6.77430    -2.90000    -2.90000    -3.15000    -4.80000    3.50000     1.77000\n\
!\n\
S       20.70290    -3.24000    -3.24000    -4.47500   -39.90000    3.50000     2.00000\n\
SM      21.30600    -3.24000    -3.24000    -4.47500   -39.90000    3.50000     1.97500\n";

}


#endif
