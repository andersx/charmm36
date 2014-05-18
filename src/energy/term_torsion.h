// term_torsion.h --- torsion angle energy term
// Copyright (C) 2009-2011 Kristoffer Enøe Johansson
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

#ifndef TERM_CHARMM36_TORSION_H
#define TERM_CHARMM36_TORSION_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "parsers/topology_parser.h"
#include "math.h"

namespace phaistos {

//! Torsion energy term
class TermCharmm36Torsion: public EnergyTermCommon<TermCharmm36Torsion, ChainFB> {

private:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36Torsion, ChainFB> EnergyTermCommon;

     std::vector<topology::DihedralAngleType9> dihedral_angles;

public:

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36Torsion(ChainFB *chain,
                        const Settings &settings=Settings(),
                        RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-torsion", settings, random_number_engine) {

          std::string filename = "/home/andersx/phaistos_dev/modules/charmm36/src/energy/charmm22_cmap/charmm22_torsion.itp";
          std::vector<topology::DihedralType9Parameter> dihedral_type_9_parameters = topology::read_dihedral_type_9_parameters(filename);

          dihedral_angles = topology::generate_dihedral_pairs(this->chain, dihedral_type_9_parameters);

     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36Torsion(const TermCharmm36Torsion &other,
                        RandomNumberEngine *random_number_engine,
                        int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            dihedral_angles(other.dihedral_angles) {}

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {


        double e_torsion = 0.0;

        for (unsigned int i = 0; i < this->dihedral_angles.size(); i++) {

            topology::DihedralAngleType9 dihedral = this->dihedral_angles[i];

            double angle = calc_dihedral((dihedral.atom1)->position,
                                         (dihedral.atom2)->position,
                                         (dihedral.atom3)->position,
                                         (dihedral.atom4)->position);

            // const double mdphi = dihedral.mult * angle - dihedral.phi0 * gmx_deg2rad;
            // const double v1 = 1.0 + cos(mdphi);
            // const double e_torsion_temp = v1 * dihedral.cp;

            const double e_torsion_temp = dihedral.cp * cos(dihedral.mult * angle - dihedral.phi0 / 180.0 * M_PI) + dihedral.cp;

            e_torsion += e_torsion_temp;

             // printf("ASC: TOR XYZ1 = %8.4f %8.4f %8.4f   XYZ2 = %8.4f %8.4f %8.4f   XYZ3 = %8.4f %8.4f %8.4f   XYZ4 = %8.4f %8.4f %8.4f   a = %15.10f   phi0 = %9.4f   cp = %8.4f   mult = %d   etor = %14.10f\n",

             //         (dihedral.atom1)->position[0], (dihedral.atom1)->position[1], (dihedral.atom1)->position[2],
             //         (dihedral.atom2)->position[0], (dihedral.atom2)->position[1], (dihedral.atom2)->position[2],
             //         (dihedral.atom3)->position[0], (dihedral.atom3)->position[1], (dihedral.atom3)->position[2],
             //         (dihedral.atom4)->position[0], (dihedral.atom4)->position[1], (dihedral.atom4)->position[2],
             //         angle *180.0 / M_PI, dihedral.phi0, dihedral.cp, dihedral.mult, e_torsion_temp);
        }

        // printf("          torsion E = %12.4f kJ/mol\n", e_torsion);
        return e_torsion / 4.184;

     }

};

}

#endif
