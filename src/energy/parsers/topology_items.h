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
#include <math.h>

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

namespace topology {


struct CmapPair {
    phaistos::Residue *residue;
    unsigned int cmap_type_index;
    int residue_index;
};


struct NonBondedParameter{

    std::string atom_type;
    unsigned int atom_number;
    double atom_mass;
    double atom_charge;
    std::string atom_ptype;
    double sigma;
    double epsilon;

};


struct NonBonded14Parameter{

    std::string atom_type1;
    std::string atom_type2;
    unsigned int pair_function;
    double sigma;
    double epsilon;

};


struct NonBondedPair{
    phaistos::Atom *atom1;
    phaistos::Atom *atom2;
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

    bool do_eef1;
    double fac_12;
    double fac_21;
    double R_vdw_1;
    double R_vdw_2;
    double lambda1;
    double lambda2;

};


struct DihedralType9Parameter {

    std::string type1;
    std::string type2;
    std::string type3;
    std::string type4;
    double phi0;
    double cp;
    unsigned int mult;

};


struct DihedralAngleType9 {

    phaistos::Atom *atom1;
    phaistos::Atom *atom2;
    phaistos::Atom *atom3;
    phaistos::Atom *atom4;
    double phi0;
    double cp;
    unsigned int mult;

};




struct BondedPairParameter {
    std::string type1;
    std::string type2;
    double kb;
    double r0;
};


struct BondedPair {
    phaistos::Atom *atom1;
    phaistos::Atom *atom2;
    double kb;
    double r0;
};


struct AngleBendParameter {
    std::string type1;
    std::string type2;
    std::string type3;
    double theta0;
    double k0;
    double r13;
    double kub;
};


struct AngleBendPair {
    phaistos::Atom *atom1;
    phaistos::Atom *atom2;
    phaistos::Atom *atom3;
    double theta0;
    double k0;
    double r13;
    double kub;
};


struct DihedralType2Parameter {

    std::string type1;
    std::string type2;
    std::string type3;
    std::string type4;
    double phi0;
    double cp;

};


struct Imptor {

    phaistos::Atom *atom1;
    phaistos::Atom *atom2;
    phaistos::Atom *atom3;
    phaistos::Atom *atom4;
    double phi0;
    double cp;

};


} // End namespace topology
#endif
