// CHARMM36 parameter parser for PHAISTOS
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

#ifndef TERM_CHARMM36_CHARMM_PARSER_H
#define TERM_CHARMM36_CHARMM_PARSER_H

#include <string>
#include <math.h>

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

namespace charmm_parser {


// Returns the atom type according to the CHARMM36 force field
// as implemented in CHARMM ... note that these are slightly
// different that the atom types defined in CHARMM36.
std::string get_charmm36_atom_type(phaistos::Atom *atom) {

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
        atom_map[CB]  = "CT2A";
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
            atom_map[CE1]  = "CPH1";
            atom_map[HE1]  = "HR3";
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
            if (res->has_atom(H3)) {
                atom_map[N]    = "NH3";
                atom_map[H1]   = "HC";
                atom_map[H2]   = "HC";
                atom_map[H3]   = "HC";
                atom_map[CA]   = "CT2";
                atom_map[HA2]  = "HB2";
                atom_map[HA3]  = "HB2";
            } else {
                atom_map[N]    = "NH2";
                atom_map[H1]   = "H";
                atom_map[H2]   = "H";
                atom_map[CA]   = "CT2";
                atom_map[HA2]  = "HB2";
                atom_map[HA3]  = "HB2";
            }
            break;
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
            if (res->has_atom(H3)) {
                atom_map[N]    = "NH3";
                atom_map[H1]   = "HC";
                atom_map[H2]   = "HC";
                atom_map[H3]   = "HC";
                atom_map[CA]   = "CT1";
                atom_map[HA]   = "HB1";
            } else {
                atom_map[N]    = "NH2";
                atom_map[H1]   = "H";
                atom_map[H2]   = "H";
                atom_map[CA]   = "CT1";
                atom_map[HA]   = "HB1";
            }
            break;
        }
    }

    if (atom->residue->terminal_status == CTERM) {
        // There should be an "if (res->has_atom(HXT))" statment here,
        // but Phaistos 1.0 does not support neutral C-term HXT atoms.
        atom_map[C]    = "CC";
        atom_map[O]    = "OC";
        atom_map[OXT]  = "OC";
    }


    if (atom_map.count(atom->atom_type) < 1)
        std::cout << "ASC: NOT found atom: " << atom << std::endl;

    return atom_map[atom->atom_type];

}

double get_charmm36_atom_charge(phaistos::Atom *atom) {

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
        atom_map[CD]  = 0.20;
        atom_map[HD2] = 0.09;
        atom_map[HD3] = 0.09;
        atom_map[NE]  = -0.70;
        atom_map[HE]  = 0.44;
        atom_map[CZ]  = 0.64;
        atom_map[NH1] = -0.80;
        atom_map[HH11]= 0.46;
        atom_map[HH12]= 0.46;
        atom_map[NH2] = -0.80;
        atom_map[HH21]= 0.46;
        atom_map[HH22]= 0.46;
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
        atom_map[HB2] = 0.09;
        atom_map[HB3] = 0.09;
        atom_map[CG]  = 0.62;
        atom_map[OD1] = -0.76;
        atom_map[OD2] = -0.76;
        atom_map[C]   = 0.51;
        atom_map[O]   = -0.51;
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
        atom_map[HG2] = 0.09;
        atom_map[HG3] = 0.09;
        atom_map[CD]  = 0.62;
        atom_map[OE1] = -0.76;
        atom_map[OE2] = -0.76;
        atom_map[C]   = 0.51;
        atom_map[O]   = -0.51;
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
        atom_map[N]    = -0.47;
        atom_map[H]    = 0.31;
        atom_map[CA]   = 0.07;
        atom_map[HA]   = 0.09;
        atom_map[CB]   = -0.05;
        atom_map[HB2]  = 0.09;
        atom_map[HB3]  = 0.09;
        atom_map[CD2]  = 0.19;
        atom_map[HD2]  = 0.13;
        atom_map[CG]   = 0.19;
        atom_map[NE2]  = -0.51;
        atom_map[HE2]  = 0.44;
        atom_map[ND1]  = -0.51;
        atom_map[HD1]  = 0.44;
        atom_map[CE1]  = 0.32;
        atom_map[HE1]  = 0.18;
        atom_map[C]    = 0.51;
        atom_map[O]    = -0.51;
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
        atom_map[CE]   = 0.21;
        atom_map[HE2]  = 0.05;
        atom_map[HE3]  = 0.05;
        atom_map[NZ]   = -0.30;
        atom_map[HZ1]  = 0.33;
        atom_map[HZ2]  = 0.33;
        atom_map[HZ3]  = 0.33;
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
        atom_map[CD1]  = 0.035;
        atom_map[HD1]  = 0.115;
        atom_map[NE1]  = -0.61;
        atom_map[HE1]  = 0.38;
        atom_map[CE2]  = 0.13;
        atom_map[CD2]  = -0.02;
        atom_map[CE3]  = -0.115;
        atom_map[HE3]  = 0.115;
        atom_map[CZ3]  = -0.115;
        atom_map[HZ3]  = 0.115;
        atom_map[CZ2]  = -0.115;
        atom_map[HZ2]  = 0.115;
        atom_map[CH2]  = -0.115;
        atom_map[HH2]  = 0.115;
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


    if (atom->residue->terminal_status == NTERM) {
        atom_map[N]    = -0.3;
        atom_map[H1]   = 0.33;
        atom_map[H2]   = 0.33;
        atom_map[H3]   = 0.33;
        atom_map[HA]   = 0.10;
        atom_map[CA]   = 0.21;
    }

    if (atom->residue->terminal_status == CTERM) {
        atom_map[C]    = 0.34;
        atom_map[O]    = -0.67;
        atom_map[OXT]  = -0.67;
    }

    if (atom_map.count(atom->atom_type) < 1)
        std::cout << "ASC: NOT found atom: " << atom << std::endl;

    return atom_map[atom->atom_type];

}

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

} // End namespace gromacs parser
#endif
