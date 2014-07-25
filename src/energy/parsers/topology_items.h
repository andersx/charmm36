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

#ifndef TERM_TOPOLOGY__H
#define TERM_TOPOLOGY__H

#include <string>

namespace topology {

//! Class to hold all information for a CMAP correction evaluation
struct CmapInteraction {
    phaistos::Residue *residue;
    unsigned int cmap_type_index;
    int residue_index;
};


//! Class to hold all parameters to evaluate the non-bonded interaction between two atoms
struct NonBondedParameter {

    std::string atom_type;
    unsigned int atom_number;
    double atom_mass;
    double atom_charge;
    std::string atom_ptype;
    double sigma;
    double epsilon;

};


//! Class to hold all parameters to evaluate the non-bonded "1-4" interaction between two atoms
struct NonBonded14Parameter {

    std::string atom_type1;
    std::string atom_type2;
    unsigned int pair_function;
    double sigma;
    double epsilon;

};


//! Class to hold all parameters (and pointer to relevant atoms) necessary to evaluate a non-bonded interaction between two atoms
struct NonBondedInteraction {

    phaistos::Atom *atom1;
    phaistos::Atom *atom2;
    double qq;
    double c6;
    double c12;
    bool is_14_interaction;
    bool do_eef1;
    double fac_12;
    double fac_21;
    double R_vdw_1;
    double R_vdw_2;
    double lambda1;
    double lambda2;

};


//! Class to hold all parameters for a torsion energy term
struct TorsionParameter {

    std::string type1;
    std::string type2;
    std::string type3;
    std::string type4;
    double phi0;
    double cp;
    unsigned int mult;

};


//! Class to hold all parameters for a torsion energy term and the four atoms involved in the interaction
struct TorsionInteraction {

    phaistos::Atom *atom1;
    phaistos::Atom *atom2;
    phaistos::Atom *atom3;
    phaistos::Atom *atom4;
    double phi0;
    double cp;
    unsigned int mult;

};


//! Class to hold parameters for a bond-stretch interaction
struct BondedPairParameter {
    std::string type1;
    std::string type2;
    double kb;
    double r0;
};


//! Class to hold parameters and pointers to two atoms in a  bond-stretch interaction
struct BondedPairInteraction {
    phaistos::Atom *atom1;
    phaistos::Atom *atom2;
    double kb;
    double r0;
};


//! Class to hold parameters for an angle-bend (and Urey-Bradley) energy term.
struct AngleBendParameter {
    std::string type1;
    std::string type2;
    std::string type3;
    double theta0;
    double k0;
    double r13;
    double kub;
};


//! Class to hold parameters for an angle-bend (and Urey-Bradley) energy term and pointers to relevant atoms
struct AngleBendInteraction {
    phaistos::Atom *atom1;
    phaistos::Atom *atom2;
    phaistos::Atom *atom3;
    double theta0;
    double k0;
    double r13;
    double kub;
};


//! Class to hold parameters for an improper torsion energy term
struct ImptorParameter {

    std::string type1;
    std::string type2;
    std::string type3;
    std::string type4;
    double phi0;
    double cp;

};


//! Class to hold parameters for an improper torsion energy term and pointers to relevant atoms
struct ImptorInteraction {

    phaistos::Atom *atom1;
    phaistos::Atom *atom2;
    phaistos::Atom *atom3;
    phaistos::Atom *atom4;
    double phi0;
    double cp;

};


} // End namespace topology
#endif
