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
//atom_type, atom_number, atom_mass, atom_charge, 
  //atom_ptype, atom_sigma, atom_epsilon;

//! Read contact map from file
//! \param chain Chain object for contact map·
std::vector<NonBondedParameter> read_nonbonded_parameters(const std::string filename) {

    // std::string filename = "/home/andersx/phaistos_dev/modules/gromacs/src/energy/charmm22_cmap/charmm22_atom_data.itp";
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
//atom_type, atom_number, atom_mass, atom_charge, 
  //atom_ptype, atom_sigma, atom_epsilon;

//! Read contact map from file
//! \param chain Chain object for contact map·
std::vector<NonBondedParameter> read_nonbonded_parameters(const std::string filename) {

    // std::string filename = "/home/andersx/phaistos_dev/modules/gromacs/src/energy/charmm22_cmap/charmm22_atom_data.itp";
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

class GromacsAtomParameter {

private:

    std::string gromacs_type;

    void set_gromacs_type(Atom *atom) {

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
            atom_map[O]    = "OT1";
            atom_map[OXT]  = "OT2";
        }


        if (atom_map.count(atom->atom_type) < 1)
            std::cout << "ASC: NOT found atom: " << atom << std::endl;
    }


public:

    std::string get_atom_type() {
        return gromacs_type;
    }

    GromacsAtomParameter(Atom *atom) {

        set_gromacs_type(atom);

    }

};




}
#endif
