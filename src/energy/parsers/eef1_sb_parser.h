// eef1_sb_parser.h -- CHARMM36/EEF1-SB atom type parameter parser for PHAISTOS
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

#ifndef EEF1_SB_PARSER_H
#define EEF1_SB_PARSER_H

#include <string>
#include <math.h>

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

namespace eef1_sb_parser {

// Returns the atom type according to the CHARMM36 force field
// as implemented in CHARMM ... note that these are slightly
// different that the atom types defined in CHARMM36 in Gromacs.

//! Returns the atom type according to the CHARMM36-EEF1-SB force field    
//! \param atom Pointer the atom for which the type is to be determined
//! \returns A string with the CHARMM36-EEF1-SB atom type
std::string get_atom_type(const phaistos::Atom *atom) {

    using namespace phaistos;
    using namespace definitions;

    Residue *res = atom->residue;
    std::map<AtomEnum, std::string> atom_map;

    switch (res->residue_type) {
    case definitions::ALA:
        atom_map[N]   = "NH1";
        atom_map[H]   = "H";
        atom_map[CA]  = "CT1";
        atom_map[HA]  = "HB1";
        atom_map[CB]  = "CT3";
        atom_map[HB1] = "HA3";
        atom_map[HB2] = "HA3";
        atom_map[HB3] = "HA3";
        atom_map[C]   = "C";
        atom_map[O]   = "O";
        break;

    case ARG:
        atom_map[N]   = "NH1";
        atom_map[H]   = "H";
        atom_map[CA]  = "CT1";
        atom_map[HA]  = "HB1";
        atom_map[CB]  = "CT2";
        atom_map[HB2] = "HA2";
        atom_map[HB3] = "HA2";
        atom_map[CG]  = "CT2";
        atom_map[HG2] = "HA2";
        atom_map[HG3] = "HA2";
        atom_map[CD]  = "CT2";
        atom_map[HD2] = "HA2";
        atom_map[HD3] = "HA2";
        atom_map[NE]  = "NC2";
        atom_map[HE]  = "HC";
        atom_map[CZ]  = "C";
        atom_map[NH1] = "NC2";
        atom_map[HH11]= "HC";
        atom_map[HH12]= "HC";
        atom_map[NH2] = "NC2";
        atom_map[HH21]= "HC";
        atom_map[HH22]= "HC";
        atom_map[C]   = "C";
        atom_map[O]   = "O";
        break;

    case ASN:
        atom_map[N]   = "NH1";
        atom_map[H]   = "H";
        atom_map[CA]  = "CT1";
        atom_map[HA]  = "HB1";
        atom_map[CB]  = "CT2";
        atom_map[HB2] = "HA2";
        atom_map[HB3] = "HA2";
        atom_map[CG]  = "CC";
        atom_map[OD1] = "O";
        atom_map[ND2] = "NH2";
        atom_map[HD21]= "H";
        atom_map[HD22]= "H";
        atom_map[C]   = "C";
        atom_map[O]   = "O";
        break;

    case ASP:
        atom_map[N]   = "NH1";
        atom_map[H]   = "H";
        atom_map[CA]  = "CT1";
        atom_map[HA]  = "HB1";
        // atom_map[CB]  = "CT2A"; // CT2A in CHARM36
        atom_map[CB]  = "CT2"; // CT2 in EEF1-SB
        atom_map[HB2] = "HA2";
        atom_map[HB3] = "HA2";
        atom_map[CG]  = "CC";
        atom_map[OD1] = "OC";
        atom_map[OD2] = "OC";
        atom_map[C]   = "C";
        atom_map[O]   = "O";

        if (res->has_atom(HD1)){
            atom_map[OD1] = "OH1";
            atom_map[OD2] = "OB";
            atom_map[HD1] = "H";
        }
        if (res->has_atom(HD2)){
            atom_map[OD1] = "OB";
            atom_map[OD2] = "OH1";
            atom_map[HD2] = "H";
        }

        break;

    case CYS:
        atom_map[N]   = "NH1";
        atom_map[H]   = "H";
        atom_map[CA]  = "CT1";
        atom_map[HA]  = "HB1";
        atom_map[CB]  = "CT2";
        atom_map[HB2] = "HA2";
        atom_map[HB3] = "HA2";
        atom_map[SG]  = "S";
        atom_map[HG] = "HS";
        atom_map[C]   = "C";
        atom_map[O]   = "O";
        break;

    case GLN:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA2";
        atom_map[HB3]  = "HA2";
        atom_map[CG]   = "CT2";
        atom_map[HG2]  = "HA2";
        atom_map[HG3]  = "HA2";
        atom_map[CD]   = "CC";
        atom_map[OE1]  = "O";
        atom_map[NE2]  = "NH2";
        atom_map[HE21] = "H";
        atom_map[HE22] = "H";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case GLU:
        atom_map[N]   = "NH1";
        atom_map[H]   = "H";
        atom_map[CA]  = "CT1";
        atom_map[HA]  = "HB1";
        atom_map[CB]  = "CT2A";
        atom_map[HB2] = "HA2";
        atom_map[HB3] = "HA2";
        atom_map[CG]  = "CT2";
        atom_map[HG2] = "HA2";
        atom_map[HG3] = "HA2";
        atom_map[CD]  = "CC";
        atom_map[OE1] = "OC";
        atom_map[OE2] = "OC";
        atom_map[C]   = "C";
        atom_map[O]   = "O";

        if (res->has_atom(HE1)){
            atom_map[OE1] = "OH1";
            atom_map[OE2] = "OB";
            atom_map[HE1] = "H";
        }
        if (res->has_atom(HE2)){
            atom_map[OE1] = "OB";
            atom_map[OE2] = "OH1";
            atom_map[HE2] = "H";
        }

        break;

    case GLY:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT2";
        atom_map[HA2]  = "HB2";
        atom_map[HA3]  = "HB2";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case HIS:
        if ( (res->has_atom(HD1)) && (res->has_atom(HE2))) {
        // HSP
            atom_map[N]    = "NH1";
            atom_map[H]    = "H";
            atom_map[CA]   = "CT1";
            atom_map[HA]   = "HB1";
            atom_map[CB]   = "CT2A";
            atom_map[HB2]  = "HA2";
            atom_map[HB3]  = "HA2";
            atom_map[CD2]  = "CPH1";
            atom_map[HD2]  = "HR1";
            atom_map[CG]   = "CPH1";
            atom_map[NE2]  = "NR3";
            atom_map[HE2]  = "H";
            atom_map[ND1]  = "NR3";
            atom_map[HD1]  = "H";
            atom_map[CE1]  = "CPH2";
            atom_map[HE1]  = "HR2";
            atom_map[C]    = "C";
            atom_map[O]    = "O";
        } else if (res->has_atom(HD1)) {
        // HSD
            atom_map[N]    = "NH1";
            atom_map[H]    = "H";
            atom_map[CA]   = "CT1";
            atom_map[HA]   = "HB1";
            atom_map[CB]   = "CT2";
            atom_map[HB2]  = "HA2";
            atom_map[HB3]  = "HA2";
            atom_map[CG]   = "CPH1";
            atom_map[ND1]  = "NR1";
            atom_map[HD1]  = "H";
            atom_map[CD2]  = "CPH1";
            atom_map[HD2]  = "HR3";
            atom_map[CE1]  = "CPH2";
            atom_map[HE1]  = "HR1";
            atom_map[NE2]  = "NR2";
            atom_map[C]    = "C";
            atom_map[O]    = "O";
        } else {
        // HSE
            atom_map[N]    = "NH1";
            atom_map[H]    = "H";
            atom_map[CA]   = "CT1";
            atom_map[HA]   = "HB1";
            atom_map[CB]   = "CT2";
            atom_map[HB2]  = "HA2";
            atom_map[HB3]  = "HA2";
            atom_map[CG]   = "CPH1";
            atom_map[CD2]  = "CPH1";
            atom_map[HD2]  = "HR3";
            atom_map[NE2]  = "NR1";
            atom_map[HE2]  = "H";
            atom_map[ND1]  = "NR2";
            atom_map[CE1]  = "CPH2";
            atom_map[HE1]  = "HR1";
            atom_map[C]    = "C";
            atom_map[O]    = "O";
        }
            break;

    case ILE:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CT1";
        atom_map[HB]   = "HA1";
        atom_map[CG2]  = "CT3";
        atom_map[HG21] = "HA3";
        atom_map[HG22] = "HA3";
        atom_map[HG23] = "HA3";
        atom_map[CG1]  = "CT2";
        atom_map[HG12] = "HA2";
        atom_map[HG13] = "HA2";
        atom_map[CD1]  = "CT3";
        atom_map[HD11] = "HA3";
        atom_map[HD12] = "HA3";
        atom_map[HD13] = "HA3";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case LEU:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA2";
        atom_map[HB3]  = "HA2";
        atom_map[CG]   = "CT1";
        atom_map[HG]   = "HA1";
        atom_map[CD1]  = "CT3";
        atom_map[HD11] = "HA3";
        atom_map[HD12] = "HA3";
        atom_map[HD13] = "HA3";
        atom_map[CD2]  = "CT3";
        atom_map[HD21] = "HA3";
        atom_map[HD22] = "HA3";
        atom_map[HD23] = "HA3";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case LYS:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA2";
        atom_map[HB3]  = "HA2";
        atom_map[CG]   = "CT2";
        atom_map[HG2]  = "HA2";
        atom_map[HG3]  = "HA2";
        atom_map[CD]   = "CT2";
        atom_map[HD2]  = "HA2";
        atom_map[HD3]  = "HA2";
        atom_map[CE]   = "CT2";
        atom_map[HE2]  = "HA2";
        atom_map[HE3]  = "HA2";
        atom_map[NZ]   = "NH3";
        atom_map[HZ1]  = "HC";
        atom_map[HZ2]  = "HC";
        atom_map[HZ3]  = "HC";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case MET:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA2";
        atom_map[HB3]  = "HA2";
        atom_map[CG]   = "CT2";
        atom_map[HG2]  = "HA2";
        atom_map[HG3]  = "HA2";
        atom_map[SD]   = "S";
        atom_map[CE]   = "CT3";
        atom_map[HE1]  = "HA3";
        atom_map[HE2]  = "HA3";
        atom_map[HE3]  = "HA3";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case PHE:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA2";
        atom_map[HB3]  = "HA2";
        atom_map[CG]   = "CA";
        atom_map[CD1]  = "CA";
        atom_map[HD1]  = "HP";
        atom_map[CE1]  = "CA";
        atom_map[HE1]  = "HP";
        atom_map[CZ]   = "CA";
        atom_map[HZ]   = "HP";
        atom_map[CD2]  = "CA";
        atom_map[HD2]  = "HP";
        atom_map[CE2]  = "CA";
        atom_map[HE2]  = "HP";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case PRO:
        atom_map[N]    = "N";
        atom_map[CD]   = "CP3";
        atom_map[HD2]  = "HA2";
        atom_map[HD3]  = "HA2";
        atom_map[CA]   = "CP1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CP2";
        atom_map[HB2]  = "HA2";
        atom_map[HB3]  = "HA2";
        atom_map[CG]   = "CP2";
        atom_map[HG2]  = "HA2";
        atom_map[HG3]  = "HA2";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case SER:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA2";
        atom_map[HB3]  = "HA2";
        atom_map[OG]   = "OH1";
        atom_map[HG]   = "H";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case THR:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CT1";
        atom_map[HB]   = "HA1";
        atom_map[OG1]  = "OH1";
        atom_map[HG1]  = "H";
        atom_map[CG2]  = "CT3";
        atom_map[HG21] = "HA3";
        atom_map[HG22] = "HA3";
        atom_map[HG23] = "HA3";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case TRP:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA2";
        atom_map[HB3]  = "HA2";
        atom_map[CG]   = "CY";
        atom_map[CD1]  = "CA";
        atom_map[HD1]  = "HP";
        atom_map[NE1]  = "NY";
        atom_map[HE1]  = "H";
        atom_map[CE2]  = "CPT";
        atom_map[CD2]  = "CPT";
        atom_map[CE3]  = "CAI";
        atom_map[HE3]  = "HP";
        atom_map[CZ3]  = "CA";
        atom_map[HZ3]  = "HP";
        atom_map[CZ2]  = "CAI";
        atom_map[HZ2]  = "HP";
        atom_map[CH2]  = "CA";
        atom_map[HH2]  = "HP";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case TYR:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA2";
        atom_map[HB3]  = "HA2";
        atom_map[CG]   = "CA";
        atom_map[CD1]  = "CA";
        atom_map[HD1]  = "HP";
        atom_map[CE1]  = "CA";
        atom_map[HE1]  = "HP";
        atom_map[CZ]   = "CA";
        atom_map[OH]   = "OH1";
        atom_map[HH]   = "H";
        atom_map[CD2]  = "CA";
        atom_map[HD2]  = "HP";
        atom_map[CE2]  = "CA";
        atom_map[HE2]  = "HP";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case VAL:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB1";
        atom_map[CB]   = "CT1";
        atom_map[HB]   = "HA1";
        atom_map[CG1]  = "CT3";
        atom_map[HG11] = "HA3";
        atom_map[HG12] = "HA3";
        atom_map[HG13] = "HA3";
        atom_map[CG2]  = "CT3";
        atom_map[HG21] = "HA3";
        atom_map[HG22] = "HA3";
        atom_map[HG23] = "HA3";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    default:
        std::cout << "ASC: Unknown residue type: " << atom->residue << std::endl;
        break;
    };

    if (atom->residue->terminal_status == NTERM) {
        switch (res->residue_type) {
        case GLY:
            atom_map[N]    = "NH3";
            atom_map[H1]   = "HC";
            atom_map[H2]   = "HC";
            atom_map[H3]   = "HC";
            atom_map[CA]   = "CT2";
            atom_map[HA2]  = "HB2";
            atom_map[HA3]  = "HB2";
            atom_map[N]    = "NH2";
            atom_map[H1]   = "H";
            atom_map[H2]   = "H";
            atom_map[CA]   = "CT2";
            atom_map[HA2]  = "HB2";
            atom_map[HA3]  = "HB2";
        case PRO:
            atom_map[N]    = "NP";
            atom_map[H1]   = "HC";
            atom_map[H2]   = "HC";
            atom_map[CD]   = "CP3";
            atom_map[HD2]  = "HA2";
            atom_map[HD3]  = "HA2";
            atom_map[CA]   = "CP1";
            atom_map[HA]   = "HB1";
            break;
        default:
            atom_map[N]    = "NH3";
            atom_map[H1]   = "HC";
            atom_map[H2]   = "HC";
            atom_map[H3]   = "HC";
            atom_map[CA]   = "CT1";
            atom_map[HA]   = "HB1";
            break;
        }
    }

    if (atom->residue->terminal_status == CTERM) {
        atom_map[C]    = "CC";
        atom_map[O]    = "OC";
        atom_map[OXT]  = "OC";
    }


    if (atom_map.count(atom->atom_type) < 1)
        std::cout << "ASC: NOT found atom: " << atom << std::endl;

    return atom_map[atom->atom_type];

}

//! Returns the atom charge according to the CHARMM36-EEF1-SB force field    
//! \param atom Pointer the atom for which the charge is to be determined
//! \returns A double containing the CHARMM36-EEF1-SB atom charge
double get_atom_charge(phaistos::Atom *atom) {

    using namespace phaistos;
    using namespace definitions;

    Residue *res = atom->residue;
    std::map<AtomEnum, double> atom_map;

    switch (res->residue_type) {
    case definitions::ALA:
        atom_map[N]   = -0.47;
        atom_map[H]   = 0.31;
        atom_map[CA]  = 0.07;
        atom_map[HA]  = 0.09;
        atom_map[CB]  = -0.27;
        atom_map[HB1] = 0.09;
        atom_map[HB2] = 0.09;
        atom_map[HB3] = 0.09;
        atom_map[C]   = 0.51;
        atom_map[O]   = -0.51;
        break;

    case ARG:
        atom_map[N]   = -0.47;
        atom_map[H]   = 0.31;
        atom_map[CA]  = 0.07;
        atom_map[HA]  = 0.09;
        atom_map[CB]  = -0.18;
        atom_map[HB2] = 0.09;
        atom_map[HB3] = 0.09;
        atom_map[CG]  = -0.18;
        atom_map[HG2] = 0.09;
        atom_map[HG3] = 0.09;
        // atom_map[CD]  = 0.20;
        atom_map[CD]  = -0.30; // EEF1-SB
        // atom_map[HD2] = 0.09;
        atom_map[HD2] = 0.05; // EEF1-SB
        // atom_map[HD3] = 0.09;
        atom_map[HD3] = 0.05; // EEF1-SB
        // atom_map[NE]  = -0.70;
        atom_map[NE]  = -0.28; // EEF1-SB
        // atom_map[HE]  = 0.44;
        atom_map[HE]  = 0.12;  // EEF1-SB
        // atom_map[CZ]  = 0.64;
        atom_map[CZ]  = -0.20; // EEF1-SB
        // atom_map[NH1] = -0.80;
        // atom_map[HH11]= 0.46;
        // atom_map[HH12]= 0.46;
        // atom_map[NH2] = -0.80;
        // atom_map[HH21]= 0.46;
        // atom_map[HH22]= 0.46;
        atom_map[NH1] = -0.121; // EEF1-SB
        atom_map[HH11]= 0.2005; // EEF1-SB
        atom_map[HH12]= 0.2005; // EEF1-SB
        atom_map[NH2] = -0.121; // EEF1-SB
        atom_map[HH21]= 0.2005; // EEF1-SB
        atom_map[HH22]= 0.2005; // EEF1-SB
        atom_map[C]   = 0.51;
        atom_map[O]   = -0.51;
        break;

    case ASN:
        atom_map[N]   =  -0.47;
        atom_map[H]   =  0.31;
        atom_map[CA]  =  0.07;
        atom_map[HA]  =  0.09;
        atom_map[CB]  =  -0.18;
        atom_map[HB2] =  0.09;
        atom_map[HB3] =  0.09;
        atom_map[CG]  =  0.55;
        atom_map[OD1] =  -0.55;
        atom_map[ND2] =  -0.62;
        atom_map[HD21]=  0.32;
        atom_map[HD22]=  0.30;
        atom_map[C]   =  0.51;
        atom_map[O]   =  -0.51;
        break;

    case ASP:
        atom_map[N]   = -0.47;
        atom_map[H]   = 0.31;
        atom_map[CA]  = 0.07;
        atom_map[HA]  = 0.09;
        atom_map[CB]  = -0.28;
        // atom_map[HB2] = 0.09; //
        // atom_map[HB3] = 0.09; //
        atom_map[HB2] = 0.14; //
        atom_map[HB3] = 0.14; //
        // atom_map[CG]  = 0.62;
        atom_map[CG]  = 1.00; //
        // atom_map[OD1] = -0.76;
        // atom_map[OD2] = -0.76;
        atom_map[OD1] = -0.50; //
        atom_map[OD2] = -0.50; //
        atom_map[C]   = 0.51;
        atom_map[O]   = -0.51;

        if (res->has_atom(HD1)){
            atom_map[CB]  = -0.21;
            atom_map[HB2] = 0.09;
            atom_map[HB3] = 0.09;
            atom_map[CG]  = 0.75;
            atom_map[OD1] = -0.61;
            atom_map[OD2] = -0.55;
            atom_map[HD1] = 0.44;
        }
        if (res->has_atom(HD2)){
            atom_map[CB]  = -0.21;
            atom_map[HB2] = 0.09;
            atom_map[HB3] = 0.09;
            atom_map[CG]  = 0.75;
            atom_map[OD1] = -0.55;
            atom_map[OD2] = -0.61;
            atom_map[HD2] = 0.44;
        }
        break;

    case CYS:
        atom_map[N]   = -0.47;
        atom_map[H]   = 0.31;
        atom_map[CA]  = 0.07;
        atom_map[HA]  = 0.09;
        atom_map[CB]  = -0.11;
        atom_map[HB2] = 0.09;
        atom_map[HB3] = 0.09;
        atom_map[SG]  = -0.23;
        atom_map[HG]  = 0.16;
        atom_map[C]   = 0.51;
        atom_map[O]   = -0.51;
        break;

    case GLN:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = -0.18;
        atom_map[HB2]  = 0.09;
        atom_map[HB3]  = 0.09;
        atom_map[CG]   = -0.18;
        atom_map[HG2]  = 0.09;
        atom_map[HG3]  = 0.09;
        atom_map[CD]   = 0.55;
        atom_map[OE1]  = -0.55;
        atom_map[NE2]  = -0.62;
        atom_map[HE21] =  0.32;
        atom_map[HE22] =  0.30;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case GLU:
        atom_map[N]   = -0.47;
        atom_map[H]   = 0.31;
        atom_map[CA]  = 0.07;
        atom_map[HA]  = 0.09;
        atom_map[CB]  = -0.18;
        atom_map[HB2] = 0.09;
        atom_map[HB3] = 0.09;
        atom_map[CG]  = -0.28;
        // atom_map[HG2] = 0.09;
        // atom_map[HG3] = 0.09;
        // atom_map[CD]  = 0.62;
        // atom_map[OE1] = -0.76;
        // atom_map[OE2] = -0.76;
        atom_map[HG2] = 0.14; //EEF1-SB specific
        atom_map[HG3] = 0.14; //EEF1-SB specific
        atom_map[CD]  = 1.00; //EEF1-SB specific
        atom_map[OE1] = -0.50;//EEF1-SB specific
        atom_map[OE2] = -0.50;//EEF1-SB specific
        atom_map[C]   = 0.51;
        atom_map[O]   = -0.51;
        if (res->has_atom(HE1)){
            atom_map[CG]  = -0.21;
            atom_map[CD]  = 0.75;
            atom_map[OE1] = -0.61;
            atom_map[OE2] = -0.55;
            atom_map[HE1] = 0.44;
        }
        if (res->has_atom(HE2)){
            atom_map[CG]  = -0.21;
            atom_map[CD]  = 0.75;
            atom_map[OE1] = -0.55;
            atom_map[OE2] = -0.61;
            atom_map[HE2] = 0.44;
        }
        break;

    case GLY:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = -0.02;
        atom_map[HA2]  = 0.09;
        atom_map[HA3]  = 0.09;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case HIS:
        if ( (res->has_atom(HD1)) && (res->has_atom(HE2))) {
        // HSP
            atom_map[N]    = -0.470;
            atom_map[H]    =  0.310;
            atom_map[CA]   =  0.070;
            atom_map[HA]   =  0.090;
            // atom_map[CB]   = -0.050;
            // atom_map[HB2]  =  0.090;
            // atom_map[HB3]  =  0.090;
            // atom_map[CD2]  =  0.190;
            // atom_map[HD2]  =  0.130;
            // atom_map[CG]   =  0.190;
            // atom_map[NE2]  = -0.510;
            // atom_map[HE2]  =  0.440;
            // atom_map[ND1]  = -0.510;
            // atom_map[HD1]  =  0.440;
            // atom_map[CE1]  =  0.320;
            // atom_map[HE1]  =  0.180;
            atom_map[CB]   = -0.100;
            atom_map[HB2]  =  0.050;
            atom_map[HB3]  =  0.050;
            atom_map[CD2]  =  0.050;
            atom_map[HD2]  =  0.000;
            atom_map[CG]   =  0.050;
            atom_map[NE2]  = -0.550;
            atom_map[HE2]  =  0.450;
            atom_map[ND1]  = -0.550;
            atom_map[HD1]  =  0.450;
            atom_map[CE1]  =  0.100;
            atom_map[HE1]  =  0.000;
            atom_map[C]    =  0.510;
            atom_map[O]    = -0.510;
        } else if (res->has_atom(HD1)) {
        // HSD
            atom_map[N]    = -0.470;
            atom_map[H]    =  0.310;
            atom_map[CA]   =  0.070;
            atom_map[HA]   =  0.090;
            atom_map[CB]   = -0.090;
            atom_map[HB2]  =  0.090;
            atom_map[HB3]  =  0.090;
            atom_map[CG]   = -0.050;
            atom_map[ND1]  = -0.360;
            atom_map[HD1]  =  0.320;
            atom_map[CD2]  =  0.220;
            atom_map[HD2]  =  0.100;
            atom_map[CE1]  =  0.250;
            atom_map[HE1]  =  0.130;
            atom_map[NE2]  = -0.700;
            atom_map[C]    =  0.510;
            atom_map[O]    = -0.510;
        } else {
        // HSE
            atom_map[N]    = -0.470;
            atom_map[H]    =  0.310;
            atom_map[CA]   =  0.070;
            atom_map[HA]   =  0.090;
            atom_map[CB]   = -0.080;
            atom_map[HB2]  =  0.090;
            atom_map[HB3]  =  0.090;
            atom_map[CG]   =  0.220;
            atom_map[CD2]  = -0.050;
            atom_map[HD2]  =  0.090;
            atom_map[NE2]  = -0.360;
            atom_map[HE2]  =  0.320;
            atom_map[ND1]  = -0.700;
            atom_map[CE1]  =  0.250;
            atom_map[HE1]  =  0.130;
            atom_map[C]    =  0.510;
            atom_map[O]    = -0.510;
        }
        break;

    case ILE:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = -0.09;
        atom_map[HB]   = 0.09;
        atom_map[CG2]  = -0.27;
        atom_map[HG21] =  0.09;
        atom_map[HG22] =  0.09;
        atom_map[HG23] =  0.09;
        atom_map[CG1]  = -0.18;
        atom_map[HG12] =  0.09;
        atom_map[HG13] =  0.09;
        atom_map[CD1]  = -0.27;
        atom_map[HD11] = 0.09;
        atom_map[HD12] = 0.09;
        atom_map[HD13] = 0.09;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case LEU:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = -0.18;
        atom_map[HB2]  = 0.09;
        atom_map[HB3]  = 0.09;
        atom_map[CG]   = -0.09;
        atom_map[HG]   = 0.09;
        atom_map[CD1]  = -0.27;
        atom_map[HD11] = 0.09;
        atom_map[HD12] = 0.09;
        atom_map[HD13] = 0.09;
        atom_map[CD2]  = -0.27;
        atom_map[HD21] = 0.09;
        atom_map[HD22] = 0.09;
        atom_map[HD23] = 0.09;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case LYS:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = -0.18;
        atom_map[HB2]  = 0.09;
        atom_map[HB3]  = 0.09;
        atom_map[CG]   = -0.18;
        atom_map[HG2]  = 0.09;
        atom_map[HG3]  = 0.09;
        atom_map[CD]   = -0.18;
        atom_map[HD2]  = 0.09;
        atom_map[HD3]  = 0.09;
        // atom_map[CE]   = 0.21;
        atom_map[CE]   = 0.00;
        // atom_map[HE2]  = 0.05;
        // atom_map[HE3]  = 0.05;
        atom_map[HE2]  = 0.00;
        atom_map[HE3]  = 0.00;
        // atom_map[NZ]   = -0.30;
        atom_map[NZ]   = -0.90;
        atom_map[HZ1]  = 0.30;
        atom_map[HZ2]  = 0.30;
        atom_map[HZ3]  = 0.30;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case MET:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = -0.18;
        atom_map[HB2]  = 0.09;
        atom_map[HB3]  = 0.09;
        atom_map[CG]   = -0.14;
        atom_map[HG2]  = 0.09;
        atom_map[HG3]  = 0.09;
        atom_map[SD]   = -0.09;
        atom_map[CE]   = -0.22;
        atom_map[HE1]  = 0.09;
        atom_map[HE2]  = 0.09;
        atom_map[HE3]  = 0.09;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case PHE:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = -0.18;
        atom_map[HB2]  = 0.09;
        atom_map[HB3]  = 0.09;
        atom_map[CG]   = 0.00;
        atom_map[CD1]  = -0.115;
        atom_map[HD1]  = 0.115;
        atom_map[CE1]  = -0.115;
        atom_map[HE1]  = 0.115;
        atom_map[CZ]   = -0.115;
        atom_map[HZ]   = 0.115;
        atom_map[CD2]  = -0.115;
        atom_map[HD2]  = 0.115;
        atom_map[CE2]  = -0.115;
        atom_map[HE2]  = 0.115;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case PRO:
        atom_map[N]    = -0.29;
        atom_map[CD]   = 0.00;
        atom_map[HD2]  = 0.09;
        atom_map[HD3]  = 0.09;
        atom_map[CA]   = 0.02;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = -0.18;
        atom_map[HB2]  = 0.09;
        atom_map[HB3]  = 0.09;
        atom_map[CG]   = -0.18;
        atom_map[HG2]  = 0.09;
        atom_map[HG3]  = 0.09;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case SER:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = 0.05;
        atom_map[HB2]  = 0.09;
        atom_map[HB3]  = 0.09;
        atom_map[OG]   = -0.66;
        atom_map[HG]   = 0.43;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case THR:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = 0.14;
        atom_map[HB]   = 0.09;
        atom_map[OG1]  = -0.66;
        atom_map[HG1]  = 0.43;
        atom_map[CG2]  = -0.27;
        atom_map[HG21] = 0.09;
        atom_map[HG22] = 0.09;
        atom_map[HG23] = 0.09;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case TRP:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = -0.18;
        atom_map[HB2]  = 0.09;
        atom_map[HB3]  = 0.09;
        atom_map[CG]   = -0.03;
        atom_map[CD1]  = -0.150;
        atom_map[HD1]  = 0.220;
        atom_map[NE1]  = -0.510;
        atom_map[HE1]  = 0.370;
        atom_map[CE2]  = 0.240;
        atom_map[CD2]  = 0.110;
        atom_map[CE3]  = -0.250;
        atom_map[HE3]  = 0.170;
        atom_map[CZ3]  = -0.200;
        atom_map[HZ3]  = 0.140;
        atom_map[CZ2]  = -0.270;
        atom_map[HZ2]  = 0.160;
        atom_map[CH2]  = -0.140;
        atom_map[HH2]  = 0.140;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case TYR:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = -0.18;
        atom_map[HB2]  = 0.09;
        atom_map[HB3]  = 0.09;
        atom_map[CG]   = 0.00;
        atom_map[CD1]  = -0.115;
        atom_map[HD1]  = 0.115;
        atom_map[CE1]  = -0.115;
        atom_map[HE1]  = 0.115;
        atom_map[CZ]   = 0.11;
        atom_map[OH]   = -0.54;
        atom_map[HH]   = 0.43;
        atom_map[CD2]  = -0.115;
        atom_map[HD2]  = 0.115;
        atom_map[CE2]  = -0.115;
        atom_map[HE2]  = 0.115;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    case VAL:
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = -0.09;
        atom_map[HB]   = 0.09;
        atom_map[CG1]  = -0.27;
        atom_map[HG11] = 0.09;
        atom_map[HG12] = 0.09;
        atom_map[HG13] = 0.09;
        atom_map[CG2]  = -0.27;
        atom_map[HG21] = 0.09;
        atom_map[HG22] = 0.09;
        atom_map[HG23] = 0.09;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
        break;

    default:
        std::cout << "ASC: Unknown residue type: " << atom->residue << std::endl;
        break;
    };


    // Neutral NTERM is currently not supported in Phaistos.
    if (atom->residue->terminal_status == NTERM) {
        switch (res->residue_type) {
        case GLY:
            atom_map[N]    = -0.9;
            atom_map[H1]   = 0.2;
            atom_map[H2]   = 0.2;
            atom_map[H3]   = 0.2;
            atom_map[HA2]  = 0.05;
            atom_map[HA3]  = 0.05;
            atom_map[CA]   = 0.20;
            break;
        case PRO:
            atom_map[N]    = -0.69;
            atom_map[H1]   = 0.2;
            atom_map[H2]   = 0.2;
            atom_map[H3]   = 0.2;
            break;
        default:
            atom_map[N]    = -0.9;
            atom_map[H1]   = 0.3;
            atom_map[H2]   = 0.3;
            atom_map[H3]   = 0.3;
            atom_map[HA]   = 0.00;
            atom_map[CA]   = 0.00;
            break;
        }

    }

    // Neutral CTERM is currently not supported in Phaistos.
    if (atom->residue->terminal_status == CTERM) {
        atom_map[C]    = 1.0;
        atom_map[O]    = -0.5;
        atom_map[OXT]  = -0.5;
    }

    if (atom_map.count(atom->atom_type) < 1) {
        std::cout << "# EEF1-SB ERROR: No charges found for atom: " << atom << std::endl;
        exit(1);
    }
    return atom_map[atom->atom_type];

}

} // End namespace gromacs parser
#endif
