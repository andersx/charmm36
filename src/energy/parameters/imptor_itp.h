// imptor_itp.h.h --- std::string'ed version of imptor.itp
// Copyright (C) 2014 Anders S. Christensen
// Copyright (c) 1991-2000, University of Groningen, The Netherlands.
// Copyright (c) 2001-2004, The GROMACS development team,
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

#ifndef TERM_CHARMM_IMPTOR_ITP_H
#define TERM_CHARMM_IMPTOR_ITP_H

namespace charmm_constants {

const static std::string imptor_itp = ";      i        j        k        l  func         phi0         kphi\n\
     HE2      HE2      CE2      CE2     2     0.000000    25.104000\n\
     HR1      NR1      NR2     CPH2     2     0.000000     4.184000\n\
     HR1      NR2      NR1     CPH2     2     0.000000     4.184000\n\
     HR3     CPH1      NR1     CPH1     2     0.000000     4.184000\n\
     HR3     CPH1      NR2     CPH1     2     0.000000     4.184000\n\
     HR3     CPH1      NR3     CPH1     2     0.000000     8.368000\n\
     HR3      NR1     CPH1     CPH1     2     0.000000     4.184000\n\
     HR3      NR2     CPH1     CPH1     2     0.000000     4.184000\n\
       N        C      CP1      CP3     2     0.000000     0.000000\n\
       C       HC       HC      NC2     2     0.000000     0.000000\n\
     NR1     CPH1     CPH2        H     2     0.000000     3.765600\n\
     NR1     CPH2     CPH1        H     2     0.000000     3.765600\n\
     NR3     CPH1     CPH2        H     2     0.000000    10.041600\n\
     NR3     CPH2     CPH1        H     2     0.000000    10.041600\n\
       O      CP1      NH2       CC     2     0.000000   376.560000\n\
       O      CT1      NH2       CC     2     0.000000   376.560000\n\
       O      CT2      NH2       CC     2     0.000000   376.560000\n\
       O      CT3      NH2       CC     2     0.000000   376.560000\n\
       O      HA1      NH2       CC     2     0.000000   376.560000\n\
       O        N      CT2       CC     2     0.000000  1004.160000\n\
       O      NH2      CP1       CC     2     0.000000   376.560000\n\
       O      NH2      CT1       CC     2     0.000000   376.560000\n\
       O      NH2      CT2       CC     2     0.000000   376.560000\n\
       O      NH2      CT3       CC     2     0.000000   376.560000\n\
       O      NH2      HA1       CC     2     0.000000   376.560000\n\
   CG2D1    CG331    NG2D1     HGA4     2     0.000000   209.200000\n\
   CG2D1    CG331    NG2P1    HGR52     2     0.000000   150.624000\n\
  CG2D1O    CG2D1    NG301     HGA4     2     0.000000   443.504000\n\
  CG2D1O    CG2D1    NG311     HGA4     2     0.000000   443.504000\n\
  CG2D1O    CG2D2    NG321     HGA4     2     0.000000   443.504000\n\
  CG2D1O    CG2D2    OG301     HGA4     2     0.000000   192.464000\n\
  CG2D1O   CG2DC1    NG301     HGA4     2     0.000000   443.504000\n\
  CG2D1O   CG2DC1    NG311     HGA4     2     0.000000   443.504000\n\
  CG2D1O   CG2DC1    OG301     HGA4     2     0.000000    83.680000\n\
  CG2D2O    CG2D1    NG301     HGA4     2     0.000000   443.504000\n\
  CG2D2O    CG2D1    NG311     HGA4     2     0.000000   443.504000\n\
  CG2D2O    CG2D2    NG321     HGA4     2     0.000000   443.504000\n\
  CG2D2O    CG2D2    OG301     HGA4     2     0.000000   192.464000\n\
  CG2D2O   CG2DC2    NG301     HGA4     2     0.000000   443.504000\n\
  CG2D2O   CG2DC2    NG311     HGA4     2     0.000000   443.504000\n\
  CG2D2O   CG2DC2    OG301     HGA4     2     0.000000    83.680000\n\
  CG2DC1   CG2R61    NG2D1     HGA4     2     0.000000   251.040000\n\
  CG2DC1   CG2DC2    NG2P1    HGR52     2     0.000000   108.784000\n\
  CG2DC2   CG2R61    NG2D1     HGA4     2     0.000000   251.040000\n\
  CG2DC2   CG2DC1    NG2P1    HGR52     2     0.000000   108.784000\n\
   CG2N1    NG321    NG321    NG2D1     2     0.000000   711.280000\n\
   CG2N1    NG2P1    NG2P1    NG2P1     2     0.000000   334.720000\n\
   CG2N1    NG2D1    NG311    NG321     2     0.000000   711.280000\n\
   CG2N2    NG2P1    NG2P1   CG2R61     2     0.000000   251.040000\n\
   CG2N2    NG2P1    NG2P1    CG331     2     0.000000   251.040000\n\
   CG2O1   CG2DC1    NG2S1    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG2DC1    NG2S2    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG2DC2    NG2S1    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG2DC2    NG2S2    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG2R61    NG2S1    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG2R61    NG2S2    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG2R62    NG2S2    OG2D1     2     0.000000  1004.160000\n\
   CG2O1    CG311    NG2S0    OG2D1     2     0.000000  1004.160000\n\
   CG2O1    CG311    NG2S1    OG2D1     2     0.000000  1004.160000\n\
   CG2O1    CG311    NG2S2    OG2D1     2     0.000000  1004.160000\n\
   CG2O1    CG311    NG311    OG2D1     2     0.000000  1004.160000\n\
   CG2O1    CG321    NG2S0    OG2D1     2     0.000000  1004.160000\n\
   CG2O1    CG321    NG2S1    OG2D1     2     0.000000  1004.160000\n\
   CG2O1    CG321    NG2S2    OG2D1     2     0.000000  1004.160000\n\
   CG2O1    CG331    NG2S0    OG2D1     2     0.000000   594.128000\n\
   CG2O1    CG331    NG2S1    OG2D1     2     0.000000  1004.160000\n\
   CG2O1    CG331    NG2S2    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG3C51    NG2S0    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG3C51    NG2S1    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG3C51    NG2S2    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG3C53    NG2S0    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG3C53    NG2S1    OG2D1     2     0.000000  1004.160000\n\
   CG2O1   CG3C53    NG2S2    OG2D1     2     0.000000  1004.160000\n\
   CG2O1    NG2S0    OG2D1    HGR52     2     0.000000   418.400000\n\
   CG2O1    NG2S2    OG2D1    HGR52     2     0.000000   552.288000\n\
   CG2O2   CG2R61    OG2D1    OG302     2     0.000000   602.496000\n\
   CG2O2    CG311    OG2D1    OG302     2     0.000000   518.816000\n\
   CG2O2    CG321    OG2D1    OG302     2     0.000000   518.816000\n\
   CG2O2    CG331    OG2D1    OG302     2     0.000000   518.816000\n\
   CG2O2   CG2R61    OG2D1    OG311     2     0.000000   443.504000\n\
   CG2O2    CG311    OG2D1    OG311     2     0.000000   543.920000\n\
   CG2O2    CG321    OG2D1    OG311     2     0.000000   543.920000\n\
   CG2O2    CG331    OG2D1    OG311     2     0.000000   543.920000\n\
   CG2O2    OG2D1    OG311    HGR52     2     0.000000   627.600000\n\
   CG2O3    OG2D2    OG2D2   CG2DC1     2     0.000000   803.328000\n\
   CG2O3    OG2D2    OG2D2   CG2DC2     2     0.000000   803.328000\n\
   CG2O3    OG2D2    OG2D2    CG2O5     2     0.000000   803.328000\n\
   CG2O3    OG2D2    OG2D2   CG2R61     2     0.000000   803.328000\n\
   CG2O3    OG2D2    OG2D2    CG301     2     0.000000   803.328000\n\
   CG2O3    OG2D2    OG2D2    CG311     2     0.000000   803.328000\n\
   CG2O3    OG2D2    OG2D2    CG314     2     0.000000   803.328000\n\
   CG2O3    OG2D2    OG2D2    CG321     2     0.000000   803.328000\n\
   CG2O3    OG2D2    OG2D2    CG331     2     0.000000   803.328000\n\
   CG2O3    OG2D2    OG2D2    HGR52     2     0.000000   560.656000\n\
   CG2O4   CG2DC1    OG2D1    HGR52     2     0.000000   117.152000\n\
   CG2O4   CG2DC2    OG2D1    HGR52     2     0.000000   117.152000\n\
   CG2O4   CG2R61    OG2D1    HGR52     2     0.000000   443.504000\n\
   CG2O4    CG321    OG2D1    HGR52     2     0.000000   418.400000\n\
   CG2O4    CG331    OG2D1    HGR52     2     0.000000   418.400000\n\
   CG2O5   CG2DC1    CG331    OG2D3     2     0.000000   736.384000\n\
   CG2O5   CG2DC2    CG331    OG2D3     2     0.000000   736.384000\n\
   CG2O5    CG2O3   CG2R61    OG2D3     2     0.000000   602.496000\n\
   CG2O5   CG2R61    CG311    OG2D3     2     0.000000   602.496000\n\
   CG2O5   CG2R61    CG321    OG2D3     2     0.000000   602.496000\n\
   CG2O5   CG2R61    CG331    OG2D3     2     0.000000   502.080000\n\
   CG2O5    CG321    CG321    OG2D3     2     0.000000   585.760000\n\
   CG2O5    CG321    CG331    OG2D3     2     0.000000   585.760000\n\
   CG2O5    CG331    CG331    OG2D3     2     0.000000   585.760000\n\
   CG2O6    NG2S2    NG2S2    OG2D1     2     0.000000   669.440000\n\
   CG2O6    OG302    OG302    OG2D1     2     0.000000  1213.360000\n\
   CG2O6    OG2D2    OG2D2    OG2D2     2     0.000000   895.376000\n\
   CG2O6    NG2S1    OG2D1    OG302     2     0.000000   518.816000\n\
   CG2O6    SG311    SG311    SG2D1     2     0.000000   669.440000\n\
  CG2R53   CG251O   NG2R53    OG2D1     2     0.000000   753.120000\n\
  CG2R53   CG252O   NG2R53    OG2D1     2     0.000000   753.120000\n\
  CG2R53   CG25C1   NG2R51    OG2D1     2     0.000000   753.120000\n\
  CG2R53   CG25C2   NG2R51    OG2D1     2     0.000000   753.120000\n\
  CG2R53   CG3C41   NG2R43    OG2D1     2     0.000000  1004.160000\n\
  CG2R53   CG3C52   NG2R53    OG2D1     2     0.000000   753.120000\n\
  CG2R53   NG2R53   NG2R53    OG2D1     2     0.000000   753.120000\n\
  CG2R53   NG2R53    OG2D1    SG311     2     0.000000   359.824000\n\
  CG2R53   NG2R53    SG2D1    SG311     2     0.000000   359.824000\n\
  CG2R63   CG2R62   NG2R61    OG2D4     2     0.000000   753.120000\n\
  CG2R63   CG2RC0   NG2R61    OG2D4     2     0.000000   753.120000\n\
  CG2R63   NG2R61   NG2R61    OG2D4     2     0.000000   753.120000\n\
  CG2R63   NG2R61   NG2R62    OG2D4     2     0.000000   753.120000\n\
  CG2R64   CG2R61   NG2R60    NG2S1     2     0.000000   158.992000\n\
  CG2R64   CG2R62   NG2R62    NG2S3     2     0.000000   502.080000\n\
  CG2R64   CG2RC0   NG2R62    NG2S3     2     0.000000   334.720000\n\
  CG2R64   NG2R61   NG2R62    NG2S3     2     0.000000   334.720000\n\
   NG2O1    OG2N1    OG2N1    CG324     2     0.000000   448.524800\n\
   NG2O1    OG2N1    OG2N1    CG334     2     0.000000   425.094400\n\
   NG2S3     HGP4     HGP4   CG2R61     2     0.000000   -20.920000\n\
   NG2S3     HGP4     HGP4   CG2R64     2     0.000000    75.312000\n\
   CC2O1    CC331    NC2D1    OC2D1     2     0.000000  1004.160000\n\
   NC2D1    CC2O1   CC3161     HCP1     2     0.000000   167.360000\n\
   NC2D1    CC2O1   CC3062     HCP1     2     0.000000   167.360000\n\
   CC2O2    CC311    OC2D2    OC2D2     2     0.000000   803.328000\n\
   CC2O2    CC301    OC2D2    OC2D2     2     0.000000   803.328000\n\
   CC2O2   CC3163    OC2D2    OC2D2     2     0.000000   803.328000\n\
   NC2D1    CC2O1    CC331     HCP1     2     0.000000   167.360000\n\
   NC2D1    CC2O1    CC311     HCP1     2     0.000000   167.360000\n\
   CC2O2    CC321    OC2D2    OC2D2     2     0.000000   803.328000\n\
   CC2O4    CC331    OC2D4     HCR1     2     0.000000   418.400000\n\
   CC2O4    CC312    OC2D4     HCR1     2     0.000000   418.400000\n\
   OC2D3    CC331    CC331    CC2O3     2     0.000000   585.760000\n\
   OC2D3    CC312    CC322    CC2O3     2     0.000000   585.760000\n\
    NN2B      CN4      CN5      HN2     2     0.000000    58.576000\n\
     NN1      CN2      HN1      HN1     2     0.000000    50.208000\n\
     CN1     NN2G     CN5G      ON1     2     0.000000   753.120000\n\
    CN1T     NN2B     NN2U      ON1     2     0.000000   920.480000\n\
     CN1     NN2U     CN3T      ON1     2     0.000000   753.120000\n\
     CN2     NN3G     NN2G      NN1     2     0.000000   334.720000\n\
     CN2     NN3A      CN5      NN1     2     0.000000   334.720000\n\
     CN2      NN3      CN3      NN1     2     0.000000   502.080000\n\
    HEL2     HEL2     CEL2     CEL2     2     0.000000    25.104000\n\
     CPB      CPA      NPH      CPA     2     0.000000   174.054400\n\
      HA      CPA      CPA      CPM     2     0.000000   246.019200\n\
     NPH      CPA      CPA       FE     2     0.000000  1149.763200\n\
     NPH      CPA      CPB      CPB     2     0.000000   339.740800\n\
     NPH      CPA      CPM      CPA     2     0.000000   153.134400\n\
     NPH      CPM      CPB      CPA     2     0.000000   273.633600\n\
     NC2        X        X        C     2     0.000000   376.560000\n\
     NC2        X        X       HC     2     0.000000   -16.736000\n\
     NH1        X        X        H     2     0.000000   167.360000\n\
     NH2        X        X        H     2     0.000000    33.472000\n\
       O        X        X        C     2     0.000000  1004.160000\n\
      OB        X        X       CD     2     0.000000   836.800000\n\
      OC        X        X       CC     2     0.000000   803.328000\n\
      CC        X        X      CT1     2     0.000000   803.328000\n\
      CC        X        X      CT2     2     0.000000   803.328000\n\
      CC        X        X      CT3     2     0.000000   803.328000\n\
     HN2        X        X      NN2     2     0.000000     8.368000\n\
     HN1        X        X      NN1     2     0.000000    33.472000\n\
     CN1        X        X      ON1     2     0.000000   753.120000\n\
    CN1T        X        X      ON1     2     0.000000   753.120000\n\
     CN1        X        X     ON1C     2     0.000000   669.440000\n\
     CN2        X        X      NN1     2     0.000000   753.120000\n\
     CN9        X        X     CN3T     2     0.000000   117.152000\n\
     OBL        X        X       CL     2     0.000000   836.800000\n\
     OCL        X        X       CL     2     0.000000   803.328000\n\
     OCL        X        X      CCL     2     0.000000   803.328000\n\
     CPB        X        X      CE1     2     0.000000   753.120000\n\
     CT2        X        X      CPB     2     0.000000   753.120000\n\
     CT3        X        X      CPB     2     0.000000   753.120000\n";

}

#endif
