// term_vdw.h --- Van der Waals interaction energy term
// Copyright (C) 2009-2014 Kristoffer Enøe Johansson, Wouter Boomsma,
// Anders S. Christensen
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

#ifndef TERM_GROMACS_PARSER_H
#define TERM_GROMACS_PARSER_H

#include <string>
#include <math.h>

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

namespace phaistos {



struct NonBondedParameter{

    std::string atom_type;
    unsigned int atom_number;
    double atom_mass;
    double atom_charge;
    std::string atom_ptype;
    double sigma;
    double epsilon;

};

std::vector<NonBondedParameter> read_nonbonded_parameters(const std::string filename) {

    std::ifstream input_stream(filename.c_str());

    if (!input_stream.is_open()) {
        std::cerr << "# Error: Cannot open UPL file " << " .\n";
        exit(EXIT_FAILURE);
    }

    std::vector<NonBondedParameter> parameters;

    while (input_stream.good()) {

        std::string line;
        std::getline(input_stream, line);

        boost::trim(line);

        if (line.size() == 0 || line[0] == ';') {
            continue;
        }

        std::vector<std::string> split_line;
        boost::split(split_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        NonBondedParameter parameter;

        parameter.atom_type   = boost::lexical_cast<std::string >(split_line[0]);
        parameter.atom_number = boost::lexical_cast<unsigned int>(split_line[1]);
        parameter.atom_mass   = boost::lexical_cast<double      >(split_line[2]);
        parameter.atom_charge = boost::lexical_cast<double      >(split_line[3]);
        parameter.atom_ptype  = boost::lexical_cast<std::string >(split_line[4]);
        parameter.sigma       = boost::lexical_cast<double      >(split_line[5]);
        parameter.epsilon     = boost::lexical_cast<double      >(split_line[6]);

        parameters.push_back(parameter);

      }

    return parameters;
}


struct NonBonded14Parameter{

    std::string atom_type1;
    std::string atom_type2;
    unsigned int pair_function;
    double sigma;
    double epsilon;

};

std::vector<NonBonded14Parameter> read_nonbonded_14_parameters(const std::string filename) {

    std::ifstream input_stream(filename.c_str());

    if (!input_stream.is_open()) {
        std::cerr << "# Error: Cannot open UPL file " << " .\n";
        exit(EXIT_FAILURE);
    }

    std::vector<NonBonded14Parameter> parameters;

    while (input_stream.good()) {

        std::string line;
        std::getline(input_stream, line);

        boost::trim(line);

        if (line.size() == 0 || line[0] == ';') {
            continue;
        }

        std::vector<std::string> split_line;
        boost::split(split_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        NonBonded14Parameter parameter;

        parameter.atom_type1    = boost::lexical_cast<std::string >(split_line[0]);
        parameter.atom_type2    = boost::lexical_cast<std::string >(split_line[1]);
        parameter.pair_function = boost::lexical_cast<unsigned int>(split_line[2]);
        parameter.sigma         = boost::lexical_cast<double      >(split_line[3]);
        parameter.epsilon       = boost::lexical_cast<double      >(split_line[4]);

        parameters.push_back(parameter);

      }

    return parameters;
}

std::string get_charmm22_atom_type(Atom *atom) {

    using namespace definitions;

    Residue *res = atom->residue;
    std::map<AtomEnum, std::string> atom_map;

    switch (res->residue_type) {
    case definitions::ALA:
        atom_map[N]   = "NH1";
        atom_map[H]   = "H";
        atom_map[CA]  = "CT1";
        atom_map[HA]  = "HB";
        atom_map[CB]  = "CT3";
        atom_map[HB1] = "HA";
        atom_map[HB2] = "HA";
        atom_map[HB3] = "HA";
        atom_map[C]   = "C";
        atom_map[O]   = "O";
        break;

    case ARG:
        atom_map[N]   = "NH1";
        atom_map[H]   = "H";
        atom_map[CA]  = "CT1";
        atom_map[HA]  = "HB";
        atom_map[CB]  = "CT2";
        atom_map[HB2] = "HA";
        atom_map[HB3] = "HA";
        atom_map[CG]  = "CT2";
        atom_map[HG2] = "HA";
        atom_map[HG3] = "HA";
        atom_map[CD]  = "CT2";
        atom_map[HD2] = "HA";
        atom_map[HD3] = "HA";
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
        atom_map[H]  = "H";
        atom_map[CA]  = "CT1";
        atom_map[HA]  = "HB";
        atom_map[CB]  = "CT2";
        atom_map[HB2] = "HA";
        atom_map[HB3] = "HA";
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
        atom_map[HA]  = "HB";
        atom_map[CB]  = "CT2";
        atom_map[HB2] = "HA";
        atom_map[HB3] = "HA";
        atom_map[CG]  = "CC";
        atom_map[OD1] = "OC";
        atom_map[OD2] = "OC";
        atom_map[C]   = "C";
        atom_map[O]   = "O";
        break;

    case CYS:
        atom_map[N]   = "NH1";
        atom_map[H]   = "H";
        atom_map[CA]  = "CT1";
        atom_map[HA]  = "HB";
        atom_map[CB]  = "CT2";
        atom_map[HB2] = "HA";
        atom_map[HB3] = "HA";
        atom_map[SG]  = "S";
        atom_map[HG1] = "HS";
        atom_map[C]   = "C";
        atom_map[O]   = "O";
        break;

    case GLN:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA";
        atom_map[HB3]  = "HA";
        atom_map[CG]   = "CT2";
        atom_map[HG2]  = "HA";
        atom_map[HG3]  = "HA";
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
        atom_map[HA]  = "HB";
        atom_map[CB]  = "CT2";
        atom_map[HB2] = "HA";
        atom_map[HB3] = "HA";
        atom_map[CG]  = "CT2";
        atom_map[HG2] = "HA";
        atom_map[HG3] = "HA";
        atom_map[CD]  = "CC";
        atom_map[OE1] = "OC";
        atom_map[OE2] = "OC";
        atom_map[C]   = "C";
        atom_map[O]   = "O";
        break;

    case GLY:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT2";
        atom_map[HA2]  = "HB";
        atom_map[HA3]  = "HB";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case HIS:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA";
        atom_map[HB3]  = "HA";
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
        break;

    case ILE:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT1";
        atom_map[HB]   = "HA";
        atom_map[CG2]  = "CT3";
        atom_map[HG21] = "HA";
        atom_map[HG22] = "HA";
        atom_map[HG23] = "HA";
        atom_map[CG1]  = "CT2";
        atom_map[HG12] = "HA";
        atom_map[HG13] = "HA";
        atom_map[CD1]  = "CT3";
        atom_map[HD11] = "HA";
        atom_map[HD12] = "HA";
        atom_map[HD13] = "HA";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case LEU:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA";
        atom_map[HB3]  = "HA";
        atom_map[CG]   = "CT1";
        atom_map[HG]   = "HA";
        atom_map[CD1]  = "CT3";
        atom_map[HD11] = "HA";
        atom_map[HD12] = "HA";
        atom_map[HD13] = "HA";
        atom_map[CD2]  = "CT3";
        atom_map[HD21] = "HA";
        atom_map[HD22] = "HA";
        atom_map[HD23] = "HA";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case LYS:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA";
        atom_map[HB3]  = "HA";
        atom_map[CG]   = "CT2";
        atom_map[HG2]  = "HA";
        atom_map[HG3]  = "HA";
        atom_map[CD]   = "CT2";
        atom_map[HD2]  = "HA";
        atom_map[HD3]  = "HA";
        atom_map[CE]   = "CT2";
        atom_map[HE2]  = "HA";
        atom_map[HE3]  = "HA";
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
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA";
        atom_map[HB3]  = "HA";
        atom_map[CG]   = "CT2";
        atom_map[HG2]  = "HA";
        atom_map[HG3]  = "HA";
        atom_map[SD]   = "S";
        atom_map[CE]   = "CT3";
        atom_map[HE1]  = "HA";
        atom_map[HE2]  = "HA";
        atom_map[HE3]  = "HA";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case PHE:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA";
        atom_map[HB3]  = "HA";
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
        atom_map[HD2]  = "HA";
        atom_map[HD3]  = "HA";
        atom_map[CA]   = "CP1";
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CP2";
        atom_map[HB2]  = "HA";
        atom_map[HB3]  = "HA";
        atom_map[CG]   = "CP2";
        atom_map[HG2]  = "HA";
        atom_map[HG3]  = "HA";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case SER:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA";
        atom_map[HB3]  = "HA";
        atom_map[OG]   = "OH1";
        atom_map[HG]   = "H";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case THR:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT1";
        atom_map[HB]   = "HA";
        atom_map[OG1]  = "OH1";
        atom_map[HG1]  = "H";
        atom_map[CG2]  = "CT3";
        atom_map[HG21] = "HA";
        atom_map[HG22] = "HA";
        atom_map[HG23] = "HA";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    case TRP:
        atom_map[N]    = "NH1";
        atom_map[H]    = "H";
        atom_map[CA]   = "CT1";
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA";
        atom_map[HB3]  = "HA";
        atom_map[CG]   = "CY";
        atom_map[CD1]  = "CA";
        atom_map[HD1]  = "HP";
        atom_map[NE1]  = "NY";
        atom_map[HE1]  = "H";
        atom_map[CE2]  = "CPT";
        atom_map[CD2]  = "CPT";
        atom_map[CE3]  = "CA";
        atom_map[HE3]  = "HP";
        atom_map[CZ3]  = "CA";
        atom_map[HZ3]  = "HP";
        atom_map[CZ2]  = "CA";
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
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT2";
        atom_map[HB2]  = "HA";
        atom_map[HB3]  = "HA";
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
        atom_map[HA]   = "HB";
        atom_map[CB]   = "CT1";
        atom_map[HB]   = "HA";
        atom_map[CG1]  = "CT3";
        atom_map[HG11] = "HA";
        atom_map[HG12] = "HA";
        atom_map[HG13] = "HA";
        atom_map[CG2]  = "CT3";
        atom_map[HG21] = "HA";
        atom_map[HG22] = "HA";
        atom_map[HG23] = "HA";
        atom_map[C]    = "C";
        atom_map[O]    = "O";
        break;

    default:
        std::cout << "ASC: Unknown residue type: " << atom->residue << std::endl;
        break;
    };

    if (atom->residue->terminal_status == NTERM) {
        atom_map[N]    = "NH3";
        atom_map[H1]   = "HC";
        atom_map[H2]   = "HC";
        atom_map[H3]   = "HC";
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

double get_charmm22_atom_charge(Atom *atom) {

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
        atom_map[HG1] = 0.16;
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
        atom_map[HA]   = 0.21;
        atom_map[CA]   = 0.10;
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


struct NonBondedPair{
    Atom *atom1;
    Atom *atom2;
    double q1;
    double q2;
    double qq;
    double sigma1;
    double sigma2;
    double epsilon1;
    double epsilon2;
    double sigma_effective;
    double epsilon_effective;
    double c6;
    double c12;
    bool is_14_interaction;
    int i1;
    int i2;
};


NonBondedParameter get_non_bonded_parameter(std::string atom_type,
        std::vector<NonBondedParameter> non_bonded_parameters) {

    NonBondedParameter parameter;

    for (unsigned int i = 0; i < non_bonded_parameters.size(); i++) {

        if (non_bonded_parameters[i].atom_type == atom_type)
            return non_bonded_parameters[i];

    }

    return parameter;

}


std::vector<NonBondedPair> generate_non_bonded_pairs(ChainFB *chain,
                    std::vector<NonBondedParameter> non_bonded_parameters,
                    std::vector<NonBonded14Parameter> non_bonded_14_parameters) {

    std::vector<NonBondedPair> non_bonded_pairs;

    int i = -1;
    for (AtomIterator<ChainFB, definitions::ALL> it1(*chain); !it1.end(); ++it1) {
        i++;
        int j = -1;

        Atom *atom1 = &*it1;

        std::string atom_type1 = get_charmm22_atom_type(atom1);
        double atom_charge1 = get_charmm22_atom_charge(atom1);

        NonBondedParameter parameter1 = get_non_bonded_parameter(atom_type1, non_bonded_parameters);

        // printf("ASC: XYZ = %7.3f  %7.3f %7.3f     q = %5.3f      ", // ix*10.0, iy*10.0, iz*10.0, charge[i]);
        //     (atom1->position)[0], (atom1->position)[1], (atom1->position)[2], atom_charge1);

        // std::cout << atom1 << std::endl;


        for (AtomIterator<ChainFB, definitions::ALL> it2(*chain); !it2.end(); ++it2) {

            j++;
            if (i > j) continue;

            Atom *atom2 = &*it2;
            int d = chain_distance<ChainFB>(atom1,atom2);
            // int test_atom_index = 295;

            std::string atom_type2 = get_charmm22_atom_type(atom2);
            double atom_charge2 = get_charmm22_atom_charge(atom2);



            NonBondedParameter parameter2 = get_non_bonded_parameter(atom_type2, non_bonded_parameters);

            if (d > 3) {
               NonBondedPair non_bonded_pair;

               double epsilon_effective = sqrt(parameter1.epsilon * parameter2.epsilon);
               double sigma_effective   = 0.5 * (parameter1.sigma + parameter2.sigma);
               // std::cout << "   " << epsilon_effective << "   " << sigma_effective << std::endl;

               non_bonded_pair.atom1 = atom1;
               non_bonded_pair.atom2 = atom2;
               non_bonded_pair.q1 = atom_charge1;
               non_bonded_pair.q2 = atom_charge2;
               non_bonded_pair.sigma1 = parameter1.sigma;
               non_bonded_pair.sigma2 = parameter2.sigma;
               non_bonded_pair.epsilon1 = parameter1.epsilon;
               non_bonded_pair.epsilon2 = parameter2.epsilon;
               non_bonded_pair.sigma_effective = sigma_effective;
               non_bonded_pair.epsilon_effective = epsilon_effective;
               non_bonded_pair.qq  = non_bonded_pair.q1 * non_bonded_pair.q2 * 138.935455;
               non_bonded_pair.c6  = 4 * epsilon_effective * std::pow(sigma_effective, 6.0);
               non_bonded_pair.c12 = 4 * epsilon_effective * std::pow(sigma_effective, 12.0);

               //  if ((std::fabs(parameter1.atom_charge - atom_charge1) > 0.001) ||
               //      (std::fabs(parameter2.atom_charge - atom_charge2) > 0.001) ) std::cout << "ERROR";
               //  std::cout << "   " << i
               //            << "   " << j
               //            << "   " << atom_type1
               //            << "   " << atom_type2
               //          //  << "   " << parameter1.sigma
               //          //  << "   " << parameter2.sigma
               //          //  << "   " << parameter1.epsilon
               //          //  << "   " << parameter2.epsilon
               //            << "   " << parameter1.atom_charge
               //            << "   " << atom_charge1
               //            << "   " << parameter2.atom_charge
               //            << "   " << atom_charge2
               //          //  << "   " << non_bonded_pair.c6
               //          //  << "   " << non_bonded_pair.c12
               //            << std::endl;


               non_bonded_pair.is_14_interaction = false;
               non_bonded_pair.i1 = i;
               non_bonded_pair.i2 = j;

               non_bonded_pairs.push_back(non_bonded_pair);

            }


            // std::cout << atom2<< i << parameter2.get_atom_type() << std::endl;

            // if ((d == 3) && ((i==test_atom_index)||(j==test_atom_index))) {
            //      std::cout << atom1 << atom2 << "   " << i << "   " << j << "   " << d << std::endl;
            // }
        }
    }

    return non_bonded_pairs;
}








}
#endif
