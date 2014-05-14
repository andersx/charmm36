// CHARMM22/CMAP parameter parser for PHAISTOS
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

NonBonded14Parameter get_non_bonded14_parameter(std::string atom_type1, std::string atom_type2,
        std::vector<NonBonded14Parameter> non_bonded14_parameters,
        std::vector<NonBondedParameter> non_bonded_parameters) {

    NonBonded14Parameter parameter;
    bool found_parameter = false;

    for (unsigned int i = 0; i < non_bonded14_parameters.size(); i++) {

        if (((non_bonded14_parameters[i].atom_type1 == atom_type1) &&
             (non_bonded14_parameters[i].atom_type2 == atom_type2)) ||
            ((non_bonded14_parameters[i].atom_type1 == atom_type2) &&
             (non_bonded14_parameters[i].atom_type2 == atom_type1))) {

            parameter = non_bonded14_parameters[i];
            found_parameter = true;
            break;

        }
    }
// struct NonBonded14Parameter{
//
//     std::string atom_type1;
//     std::string atom_type2;
//     unsigned int pair_function;
//     double sigma;
//     double epsilon;
//
// };


    if (!found_parameter) {

        NonBondedParameter parameter1 = get_non_bonded_parameter(atom_type1, non_bonded_parameters);
        NonBondedParameter parameter2 = get_non_bonded_parameter(atom_type2, non_bonded_parameters);

        parameter.atom_type1 = atom_type1;
        parameter.atom_type2 = atom_type1;

        const double sigma1 = parameter1.sigma;
        const double sigma2 = parameter2.sigma;
        const double epsilon1 = parameter1.epsilon;
        const double epsilon2 = parameter2.epsilon;

        const double epsilon_effective = sqrt(parameter1.epsilon * parameter2.epsilon);
        const double sigma_effective   = 0.5 * (parameter1.sigma + parameter2.sigma);

        parameter.sigma = sigma_effective;
        parameter.epsilon = epsilon_effective;

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

               non_bonded_pair.is_14_interaction = false;
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
//               } else if (d == 3) {
//               NonBondedPair non_bonded_pair;
//
//               double epsilon_effective = sqrt(parameter1.epsilon * parameter2.epsilon);
//               double sigma_effective   = 0.5 * (parameter1.sigma + parameter2.sigma);
//               // std::cout << "   " << epsilon_effective << "   " << sigma_effective << std::endl;
//
//               non_bonded_pair.atom1 = atom1;
//               non_bonded_pair.atom2 = atom2;
//               non_bonded_pair.q1 = atom_charge1;
//               non_bonded_pair.q2 = atom_charge2;
//               non_bonded_pair.sigma1 = parameter1.sigma;
//               non_bonded_pair.sigma2 = parameter2.sigma;
//               non_bonded_pair.epsilon1 = parameter1.epsilon;
//               non_bonded_pair.epsilon2 = parameter2.epsilon;
//               non_bonded_pair.sigma_effective = sigma_effective;
//               non_bonded_pair.epsilon_effective = epsilon_effective;
//               non_bonded_pair.qq  = non_bonded_pair.q1 * non_bonded_pair.q2 * 138.935455;
//               non_bonded_pair.c6  = 4 * epsilon_effective * std::pow(sigma_effective, 6.0);
//               non_bonded_pair.c12 = 4 * epsilon_effective * std::pow(sigma_effective, 12.0);
//
//               non_bonded_pair.is_14_interaction = false;
//               //  if ((std::fabs(parameter1.atom_charge - atom_charge1) > 0.001) ||
//               //      (std::fabs(parameter2.atom_charge - atom_charge2) > 0.001) ) std::cout << "ERROR";
//               //  std::cout << "   " << i
//               //            << "   " << j
//               //            << "   " << atom_type1
//               //            << "   " << atom_type2
//               //          //  << "   " << parameter1.sigma
//               //          //  << "   " << parameter2.sigma
//               //          //  << "   " << parameter1.epsilon
//               //          //  << "   " << parameter2.epsilon
//               //            << "   " << parameter1.atom_charge
//               //            << "   " << atom_charge1
//               //            << "   " << parameter2.atom_charge
//               //            << "   " << atom_charge2
//               //          //  << "   " << non_bonded_pair.c6
//               //          //  << "   " << non_bonded_pair.c12
//               //            << std::endl;
//
//
//               non_bonded_pair.is_14_interaction = true;
//               non_bonded_pair.i1 = i;
//               non_bonded_pair.i2 = j;
//
//               non_bonded_pairs.push_back(non_bonded_pair);
//               }
             } else if (d == 3) {


                NonBonded14Parameter parameter14 = get_non_bonded14_parameter(atom_type1, atom_type2,
                                                                            non_bonded_14_parameters,
                                                                            non_bonded_parameters);

                NonBondedParameter parameter2 = get_non_bonded_parameter(atom_type2, non_bonded_parameters);
                NonBondedPair non_bonded_pair;

                non_bonded_pair.atom1 = atom1;
                non_bonded_pair.atom2 = atom2;
                non_bonded_pair.q1 = atom_charge1;
                non_bonded_pair.q2 = atom_charge2;
                non_bonded_pair.sigma1 = parameter1.sigma;
                non_bonded_pair.sigma2 = parameter2.sigma;
                non_bonded_pair.epsilon1 = parameter1.epsilon;
                non_bonded_pair.epsilon2 = parameter2.epsilon;
                non_bonded_pair.sigma_effective = parameter14.sigma;
                non_bonded_pair.epsilon_effective = parameter14.epsilon;
                non_bonded_pair.qq  = non_bonded_pair.q1 * non_bonded_pair.q2 * 138.935455;
                non_bonded_pair.c6  = 4 * parameter14.epsilon * std::pow(parameter14.sigma, 6.0);
                non_bonded_pair.c12 = 4 * parameter14.epsilon * std::pow(parameter14.sigma, 12.0);

                non_bonded_pair.is_14_interaction = true;

                non_bonded_pairs.push_back(non_bonded_pair);
         }


// struct NonBonded14Parameter{
//
//     std::string atom_type1;
//     std::string atom_type2;
//     unsigned int pair_function;
//     double sigma;
//     double epsilon;
//
// };
//
            // std::cout << atom2<< i << parameter2.get_atom_type() << std::endl;

            // if ((d == 3) && ((i==test_atom_index)||(j==test_atom_index))) {
            //      std::cout << atom1 << atom2 << "   " << i << "   " << j << "   " << d << std::endl;
            // }
        }
    }

    return non_bonded_pairs;
}


struct DihedralType9Parameter {

    std::string type1;
    std::string type2;
    std::string type3;
    std::string type4;
    double phi0;
    double cp;
    unsigned int mult;

};


std::vector<DihedralType9Parameter> read_dihedral_type_9_parameters(const std::string filename) {

    std::ifstream input_stream(filename.c_str());

    if (!input_stream.is_open()) {
        std::cerr << "# Error: Cannot open itp file " << " .\n";
        exit(EXIT_FAILURE);
    }

    std::vector<DihedralType9Parameter> parameters;

    while (input_stream.good()) {

        std::string line;
        std::getline(input_stream, line);

        boost::trim(line);

        // std::cout << line << std::endl;
        if (line.size() == 0 || line[0] == ';') {
            continue;
        }

        std::vector<std::string> split_line;
        boost::split(split_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        DihedralType9Parameter parameter;

        parameter.type1 = boost::lexical_cast<std::string >(split_line[0]);
        parameter.type2 = boost::lexical_cast<std::string >(split_line[1]);
        parameter.type3 = boost::lexical_cast<std::string >(split_line[2]);
        parameter.type4 = boost::lexical_cast<std::string >(split_line[3]);
        parameter.phi0  = boost::lexical_cast<double>(split_line[5]);
        parameter.cp    = boost::lexical_cast<double>(split_line[6]);
        parameter.mult  = boost::lexical_cast<unsigned int >(split_line[7]);

        parameters.push_back(parameter);

      }

    return parameters;
}


struct DihedralAngleType9 {

    Atom *atom1;
    Atom *atom2;
    Atom *atom3;
    Atom *atom4;
    double phi0;
    double cp;
    unsigned int mult;

};


std::vector<DihedralAngleType9> generate_non_bonded_pairs(ChainFB *chain,
                    std::vector<DihedralType9Parameter> dihedral_angle_type_9_parameters) {

    std::vector<DihedralAngleType9> dihedral_angle_type_9s;

    // all torsions
    for (AtomIterator<ChainFB, definitions::ALL> it2(*chain); !it2.end(); ++it2) {
        Atom *atom2 = &*it2;
        for (AtomIterator<ChainFB, definitions::ALL> it3(*chain); !it3.end(); ++it3) {
            Atom *atom3 = &*it3;

            if ((atom2->residue->index >
                 atom3->residue->index ) ||
                ((atom2->residue->index == atom3->residue->index) &&
                 (atom2->index > atom3->index)))
                continue;

            if(chain_distance<ChainFB>(atom2,atom3) != 1)
                continue;
            // this->counter++;
            // std::cout << counter << atom2 << atom3 << std::endl;
            for (CovalentBondIterator<ChainFB> it1(atom2, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
                !it1.end(); ++it1) {
                Atom *atom1 = &*it1;
                for (CovalentBondIterator<ChainFB> it4(atom3, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
                    !it4.end(); ++it4) {
                    Atom *atom4 = &*it4;

                    if ((atom2 != atom4) && (atom1 != atom3)) {

                        bool found_this_parameter = false;

                        // std::cout << atom1 << atom2 << atom3 << atom4 << std::endl;

                        std::string type1 = get_charmm22_atom_type(atom1);
                        std::string type2 = get_charmm22_atom_type(atom2);
                        std::string type3 = get_charmm22_atom_type(atom3);
                        std::string type4 = get_charmm22_atom_type(atom4);

                        for (unsigned int i = 0; i < dihedral_angle_type_9_parameters.size(); i++) {

                            DihedralType9Parameter p = dihedral_angle_type_9_parameters[i];
                            // std::cout << p.type1 << type1 << std::endl;

                            if ((p.type1 == type1) && (p.type2 == type2) && (p.type3 == type3)&& (p.type4 == type4)) {
                                DihedralAngleType9 dihedral_angle_type_9;
                                dihedral_angle_type_9.atom1 = atom1;
                                dihedral_angle_type_9.atom2 = atom2;
                                dihedral_angle_type_9.atom3 = atom3;
                                dihedral_angle_type_9.atom4 = atom4;
                                dihedral_angle_type_9.phi0  = p.phi0;
                                dihedral_angle_type_9.cp    = p.cp;
                                dihedral_angle_type_9.mult  = p.mult;
                                found_this_parameter = true;

                                dihedral_angle_type_9s.push_back(dihedral_angle_type_9);
                            }
                            if ((p.type1 == type4) && (p.type2 == type3) && (p.type3 == type2)&& (p.type4 == type1)) {
                                DihedralAngleType9 dihedral_angle_type_9;
                                dihedral_angle_type_9.atom1 = atom4;
                                dihedral_angle_type_9.atom2 = atom3;
                                dihedral_angle_type_9.atom3 = atom2;
                                dihedral_angle_type_9.atom4 = atom1;
                                dihedral_angle_type_9.phi0  = p.phi0;
                                dihedral_angle_type_9.cp    = p.cp;
                                dihedral_angle_type_9.mult  = p.mult;
                                found_this_parameter = true;
                                dihedral_angle_type_9s.push_back(dihedral_angle_type_9);
                            }
                        }
                        if (!found_this_parameter) {
                            for (unsigned int i = 0; i < dihedral_angle_type_9_parameters.size(); i++) {

                                DihedralType9Parameter p = dihedral_angle_type_9_parameters[i];
                                // std::cout << p.type1 << type1 << std::endl;

                                if ((p.type1 == "X") && (p.type2 == type2) && (p.type3 == type3)&& (p.type4 == "X")) {
                                    DihedralAngleType9 dihedral_angle_type_9;
                                    dihedral_angle_type_9.atom1 = atom1;
                                    dihedral_angle_type_9.atom2 = atom2;
                                    dihedral_angle_type_9.atom3 = atom3;
                                    dihedral_angle_type_9.atom4 = atom4;
                                    dihedral_angle_type_9.phi0  = p.phi0;
                                    dihedral_angle_type_9.cp    = p.cp;
                                    dihedral_angle_type_9.mult  = p.mult;
                                    found_this_parameter = true;
                                    dihedral_angle_type_9s.push_back(dihedral_angle_type_9);
                                    break;
                                }
                                if ((p.type1 == "X") && (p.type2 == type3) && (p.type3 == type2)&& (p.type4 == "X")) {
                                    DihedralAngleType9 dihedral_angle_type_9;
                                    dihedral_angle_type_9.atom1 = atom4;
                                    dihedral_angle_type_9.atom2 = atom3;
                                    dihedral_angle_type_9.atom3 = atom2;
                                    dihedral_angle_type_9.atom4 = atom1;
                                    dihedral_angle_type_9.phi0  = p.phi0;
                                    dihedral_angle_type_9.cp    = p.cp;
                                    dihedral_angle_type_9.mult  = p.mult;
                                    found_this_parameter = true;
                                    dihedral_angle_type_9s.push_back(dihedral_angle_type_9);
                                    break;
                                }
                            }
                            if (!found_this_parameter) {
                                std::cout << "ASC: TORERR COULD NOT FIND DIHEDRAL PARAM" << atom1 << atom2 << atom3 << atom4 << std::endl;
                                std::cout << "ASC: TORERR COULD NOT FIND DIHEDRAL PARAM   " << type1 << "  " << type2 << "  " << type3 << "  " << type4 << "  " << std::endl;
                            } else {

                                // std::cout << "ASC: TOR FOUND DIHEDRAL PARAM" << atom1 << atom2 << atom3 << atom4 << std::endl;
                            }

                            // std::cout << "ASC: TORERR COULD NOT FIND DIHEDRAL PARAM" << atom1 << atom2 << atom3 << atom4 << std::endl;
                            // std::cout << "ASC: TORERR COULD NOT FIND DIHEDRAL PARAM   " << type1 << "  " << type2 << "  " << type3 << "  " << type4 << "  " << std::endl;
                        } else {

                            // std::cout << "ASC: TOR FOUND DIHEDRAL PARAM" << atom1 << atom2 << atom3 << atom4 << std::endl;
                        }
                    }
                }
            }
        }
    }
    return dihedral_angle_type_9s;
}


struct BondedPairParameter {
    std::string type1;
    std::string type2;
    double kb;
    double r0;
};


std::vector<BondedPairParameter> read_bonded_pair_parameters(const std::string filename) {

    std::ifstream input_stream(filename.c_str());

    if (!input_stream.is_open()) {
        std::cerr << "# Error: Cannot open itp file " << " .\n";
        exit(EXIT_FAILURE);
    }

    std::vector<BondedPairParameter> parameters;

    while (input_stream.good()) {

        std::string line;
        std::getline(input_stream, line);

        boost::trim(line);

        // std::cout << line << std::endl;
        if (line.size() == 0 || line[0] == ';') {
            continue;
        }

        std::vector<std::string> split_line;
        boost::split(split_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        BondedPairParameter parameter;

        parameter.type1 = boost::lexical_cast<std::string >(split_line[0]);
        parameter.type2 = boost::lexical_cast<std::string >(split_line[1]);
        parameter.r0    = boost::lexical_cast<double>(split_line[3]);
        parameter.kb    = boost::lexical_cast<double>(split_line[4]);

        parameters.push_back(parameter);

      }

    return parameters;
}


struct BondedPair {
    Atom *atom1;
    Atom *atom2;
    double kb;
    double r0;
};



std::vector<BondedPair> generate_bonded_pairs(ChainFB *chain,
                    std::vector<BondedPairParameter> bonded_pair_parameters) {

    std::vector<BondedPair> bonded_pairs;


    for (AtomIterator<ChainFB,definitions::ALL> it1(*(chain)); !it1.end(); ++it1) {
        Atom *atom1 = &*it1;
        std::string type1 = get_charmm22_atom_type(atom1);

        for (CovalentBondIterator<ChainFB> it2(atom1, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
            !it2.end(); ++it2) {

            Atom *atom2 = &*it2;
            if (atom1->residue->index < atom2->residue->index ||
                (atom1->residue->index == atom2->residue->index && atom1->index < atom2->index)) {

                std::string type2 = get_charmm22_atom_type(atom2);

                for (unsigned int i = 0; i < bonded_pair_parameters.size(); i++) {

                    BondedPairParameter parameter = bonded_pair_parameters[i];

                    if (((parameter.type1 == type1) && (parameter.type2 == type2)) ||
                        ((parameter.type1 == type2) && (parameter.type2 == type1))) {

                        BondedPair pair;
                        pair.atom1 = atom1;
                        pair.atom2 = atom2;
                        pair.kb = parameter.kb;
                        pair.r0 = parameter.r0;

                        bonded_pairs.push_back(pair);

                        break;
                    }
                }
            }
        }
    }

    return bonded_pairs;

}



struct AngleBendParameter {
    std::string type1;
    std::string type2;
    std::string type3;
    double theta0;
    double k0;
    double r13;
    double kub;
};


std::vector<AngleBendParameter> read_angle_bend_parameters(const std::string filename) {

    std::ifstream input_stream(filename.c_str());

    if (!input_stream.is_open()) {
        std::cerr << "# Error: Cannot open itp file " << " .\n";
        exit(EXIT_FAILURE);
    }

    std::vector<AngleBendParameter> parameters;

    while (input_stream.good()) {

        std::string line;
        std::getline(input_stream, line);

        boost::trim(line);

        // std::cout << line << std::endl;
        if (line.size() == 0 || line[0] == ';') {
            continue;
        }

        std::vector<std::string> split_line;
        boost::split(split_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        AngleBendParameter parameter;

        parameter.type1  = boost::lexical_cast<std::string >(split_line[0]);
        parameter.type2  = boost::lexical_cast<std::string >(split_line[1]);
        parameter.type3  = boost::lexical_cast<std::string >(split_line[2]);
        parameter.theta0 = boost::lexical_cast<double>(split_line[4]);
        parameter.k0     = boost::lexical_cast<double>(split_line[5]);
        parameter.r13    = boost::lexical_cast<double>(split_line[6]);
        parameter.kub    = boost::lexical_cast<double>(split_line[7]);

        parameters.push_back(parameter);

      }

    return parameters;
}



struct AngleBendPair {
    Atom *atom1;
    Atom *atom2;
    Atom *atom3;
    double theta0;
    double k0;
    double r13;
    double kub;
};


std::vector<AngleBendPair> generate_angle_bend_pairs(ChainFB *chain,
                    std::vector<AngleBendParameter> angle_bend_parameters) {

    std::vector<AngleBendPair> angle_bend_pairs;

    for (AtomIterator<ChainFB,definitions::ALL> it1(*(chain)); !it1.end(); ++it1) {
        Atom *atom2 = &*it1;

        std::string type2 = get_charmm22_atom_type(atom2);
        for (CovalentBondIterator<ChainFB> it2(atom2, CovalentBondIterator<ChainFB>::DEPTH_1_ONLY);
            !it2.end(); ++it2) {

            Atom *atom1 = &*it2;
            std::string type1 = get_charmm22_atom_type(atom1);

            CovalentBondIterator<ChainFB> it3(it2);
            // Fancy way to discard it3 = it1
            ++it3;
            for (; !it3.end(); ++it3) {
                Atom *atom3 = &*it3;

                std::string type3 = get_charmm22_atom_type(atom3);

                bool found_this_parameter = false;

                for (unsigned int i = 0; i < angle_bend_parameters.size(); i++) {

                    AngleBendParameter parameter = angle_bend_parameters[i];

                    if (((parameter.type1 == type1) && (parameter.type2 == type2) && (parameter.type3 == type3)) ||
                        ((parameter.type1 == type3) && (parameter.type2 == type2) && (parameter.type3 == type1))) {

                        AngleBendPair pair;
                        pair.atom1 = atom1;
                        pair.atom2 = atom2;
                        pair.atom3 = atom3;
                        pair.theta0 = parameter.theta0;
                        pair.k0     = parameter.k0;
                        pair.r13    = parameter.r13;
                        pair.kub    = parameter.kub;

                        found_this_parameter = true;
                        angle_bend_pairs.push_back(pair);

                        break;
                    }
                }
                if (!found_this_parameter) {
                    std::cout << "ASC: COULD NOT FIND angle bend PARAM" << atom1 << atom2 << atom3 << std::endl;
                    std::cout << "ASC: COULD NOT FIND angle-bend PARAM   " << type1 << "  " << type2 << "  " << type3 << std::endl;
                }
            }
        }
    }

    return angle_bend_pairs;

}



struct DihedralType2Parameter {

    std::string type1;
    std::string type2;
    std::string type3;
    std::string type4;
    double phi0;
    double cp;

};


std::vector<DihedralType2Parameter> read_dihedral_type_2_parameters(const std::string filename) {

    std::ifstream input_stream(filename.c_str());

    if (!input_stream.is_open()) {
        std::cerr << "# Error: Cannot open itp file " << " .\n";
        exit(EXIT_FAILURE);
    }

    std::vector<DihedralType2Parameter> parameters;

    while (input_stream.good()) {

        std::string line;
        std::getline(input_stream, line);

        boost::trim(line);

        // std::cout << line << std::endl;
        if (line.size() == 0 || line[0] == ';') {
            continue;
        }

        std::vector<std::string> split_line;
        boost::split(split_line, line, boost::is_any_of(" \t"), boost::token_compress_on);

        DihedralType2Parameter parameter;

        parameter.type1 = boost::lexical_cast<std::string >(split_line[0]);
        parameter.type2 = boost::lexical_cast<std::string >(split_line[1]);
        parameter.type3 = boost::lexical_cast<std::string >(split_line[2]);
        parameter.type4 = boost::lexical_cast<std::string >(split_line[3]);
        parameter.phi0  = boost::lexical_cast<double>(split_line[5]);
        parameter.cp    = boost::lexical_cast<double>(split_line[6]);

        parameters.push_back(parameter);

      }

    return parameters;
}


struct Imptor {

    Atom *atom1;
    Atom *atom2;
    Atom *atom3;
    Atom *atom4;
    double phi0;
    double cp;

};


Imptor atoms_to_imptor(std::vector<Atom*> atoms,
    std::vector<DihedralType2Parameter> imptor_parameters) {

    Imptor imptor;

    imptor.atom1 = atoms[0];
    imptor.atom2 = atoms[1];
    imptor.atom3 = atoms[2];
    imptor.atom4 = atoms[3];

    std::string type1 = get_charmm22_atom_type(atoms[0]);
    std::string type2 = get_charmm22_atom_type(atoms[1]);
    std::string type3 = get_charmm22_atom_type(atoms[2]);
    std::string type4 = get_charmm22_atom_type(atoms[3]);

    bool found_parameter = false;

    for (unsigned int i = 0; i < imptor_parameters.size(); i++) {

        DihedralType2Parameter p = imptor_parameters[i];

        if (((p.type1 == type1) && (p.type2 == type2) && (p.type3 == type3) && (p.type4 == type4)) ||
            ((p.type1 == type4) && (p.type2 == type3) && (p.type3 == type2) && (p.type4 == type1)) ||
            ((p.type1 == type1) && (p.type2 == "X"  ) && (p.type3 == "X"  ) && (p.type4 == type4)) ||
            ((p.type1 == type4) && (p.type2 == "X"  ) && (p.type3 == "X"  ) && (p.type4 == type1))) {

            found_parameter = true;

            imptor.phi0 = p.phi0;
            imptor.cp = p.cp;
        }
    }

    if (!found_parameter)
        std::cout << "ERROR: DIDN'T FIND PARAMETER" << std::endl;

    return imptor;


}


std::vector<Imptor> generate_imptors(ChainFB *chain,
    std::vector<DihedralType2Parameter> imptor_parameters) {

    using namespace definitions;
    using namespace vector_utils;

    std::vector<Imptor> imptors;

    for (ResidueIterator<ChainFB> res(*(chain)); !(res).end(); ++res) {


        std::vector<std::vector<AtomEnum> > enum_pairs;


        switch (res->residue_type) {
        case ALA:
            break;

        case ARG:
            enum_pairs.push_back(make_vector(CZ,      NH1,     NH2,     NE));
            break;

        case ASN:
            enum_pairs.push_back(make_vector(CG,      ND2,     CB,      OD1));
            enum_pairs.push_back(make_vector(CG,      CB,      ND2,     OD1));
            enum_pairs.push_back(make_vector(ND2,     CG,      HD21,    HD22));
            enum_pairs.push_back(make_vector(ND2,    CG,      HD22,    HD21));
            break;

        case ASP:
            enum_pairs.push_back(make_vector(CG,      CB,      OD2,     OD1));
            break;

        case CYS:
            break;

        case GLN:
	        enum_pairs.push_back(make_vector(CD,	NE2,	CG,	OE1));
        	enum_pairs.push_back(make_vector(CD,	CG,	NE2,	OE1));
        	enum_pairs.push_back(make_vector(NE2,	CD,	HE21,	HE22));
        	enum_pairs.push_back(make_vector(NE2,	CD,	HE22,	HE21));
            break;

        case GLU:
         	enum_pairs.push_back(make_vector(CD,	CG,	OE2,	OE1));
            break;

        case GLY:
            break;

        case HIS:
            if (res->has_atom(HD1) && res->has_atom(HE2)) {
            	enum_pairs.push_back(make_vector(ND1, 	CG	, CE1,	HD1));
            	enum_pairs.push_back(make_vector(ND1, 	CE1	, CG , HD1));
            	enum_pairs.push_back(make_vector(NE2, 	CD2	, CE1,	HE2));
            	enum_pairs.push_back(make_vector(NE2, 	CE1	, CD2,	HE2));
            } else if (!(res->has_atom(HD1)) && res->has_atom(HE2)) {
            	enum_pairs.push_back(make_vector(NE2, 	CD2	, CE1,	HE2));
            	enum_pairs.push_back(make_vector(CD2, 	CG	, NE2,	HD2));
            	enum_pairs.push_back(make_vector(CE1, 	ND1	, NE2,	HE1));
            	enum_pairs.push_back(make_vector(NE2, 	CE1	, CD2,	HE2));
            	enum_pairs.push_back(make_vector(CD2, 	NE2	, CG , HD2));
        	    enum_pairs.push_back(make_vector(CE1, 	NE2	, ND1,	HE1));
            } else if (res->has_atom(HD1) && !(res->has_atom(HE2))) {
            	enum_pairs.push_back(make_vector(ND1, 	CG	, CE1,	HD1));
            	enum_pairs.push_back(make_vector(CD2, 	CG	, NE2,	HD2));
            	enum_pairs.push_back(make_vector(CE1, 	ND1	, NE2,	HE1));
            	enum_pairs.push_back(make_vector(ND1, 	CE1	, CG , HD1));
        	    enum_pairs.push_back(make_vector(CD2, 	NE2	, CG , HD2));
            	enum_pairs.push_back(make_vector(CE1, 	NE2	, ND1,	HE1));
            } else {
                std::cout << "Unknown protonations state on " << res << std::endl;
            }
            break;

        case ILE:
            break;

        case LEU:
            break;

        case LYS:
            break;

        case MET:
            break;

        case PHE:
            break;

        case PRO:
            break;

        case SER:
            break;

        case THR:
            break;

        case TRP:
            break;

        case TYR:
            break;

        case VAL:
            break;

        default:
            std::cout << "ASC: Unknown residue type: " << res << std::endl;
            break;
        };


        for (unsigned int i = 0; i < enum_pairs.size(); i++) {

            Imptor sc_imptor = atoms_to_imptor(make_vector((*res)[enum_pairs[i][0]],
                                                           (*res)[enum_pairs[i][1]],
                                                           (*res)[enum_pairs[i][2]],
                                                           (*res)[enum_pairs[i][3]]),
                                               imptor_parameters);

            imptors.push_back(sc_imptor);


        }


        Residue *previous_residue = res->get_neighbour(-1);
        Residue *next_residue     = res->get_neighbour(+1);

        // std::cout << *res << std::endl;

         //N   -C  CA  HN
         if (!(res->terminal_status == NTERM)) {

             AtomEnum amide_atom = H;

             if (res->residue_type == PRO)
                amide_atom = CD;

             Imptor bb_imptor = atoms_to_imptor(make_vector((*res)[N],
                                                            (*previous_residue)[C],
                                                            (*res)[CA],
                                                            (*res)[amide_atom]),
                                                imptor_parameters);

             imptors.push_back(bb_imptor);
         }

          //C   CA  +N  O
         if (!(res->terminal_status == CTERM)) {

             Imptor bb_imptor = atoms_to_imptor(make_vector((*res)[C],
                                                            (*res)[CA],
                                                            (*next_residue)[N],
                                                            (*res)[O]),
                                                imptor_parameters);

             imptors.push_back(bb_imptor);
         } else if (res->terminal_status == CTERM) {

             Imptor bb_imptor = atoms_to_imptor(make_vector((*res)[C],
                                                            (*res)[CA],
                                                            (*res)[OXT],
                                                            (*res)[O]),
                                                imptor_parameters);

             imptors.push_back(bb_imptor);
         }










    }

    return imptors;

}

} // End namespace phaistos
#endif
