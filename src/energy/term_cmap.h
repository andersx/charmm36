// term_cmap.h --- CMAP torsion angle energy term
// Copyright (C) 2014 Anders S. Christensen
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

#ifndef TERM_CHARMM36_CMAP_H
#define TERM_CHARMM36_CMAP_H

#include <boost/type_traits/is_base_of.hpp>
#include "energy/energy_term.h"
#include "protein/iterators/pair_iterator_chaintree.h"

#include "protein/iterators/pair_iterator_chaintree.h"
#include "protein/chain_fb.h"
#include "protein/definitions.h"

#include "parsers/topology_parser.h"
#include "kernels/kernel_cmap.h"
#include "term_cmap_tables.h"

namespace phaistos {

//! CMAP energy term
class TermCharmm36Cmap: public EnergyTermCommon<TermCharmm36Cmap, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36Cmap, ChainFB> EnergyTermCommon;

public:


     //! Table which contains the CMAP correction tables
     //! in Gromacs' 
     std::vector<std::vector<double> > cmapdata;

     //! Vector containing all terms in the CMAP correction.
     std::vector<topology::CmapPair> cmap_pairs;

     //! Use same settings as base class
     typedef EnergyTerm<ChainFB>::SettingsClassicEnergy Settings;

     //! Constructor.
     //! \param chain Molecule chain
     //! \param settings Local Settings object
     //! \param random_number_engine Object from which random number generators can be created.
     TermCharmm36Cmap(ChainFB *chain,
                        const Settings &settings=Settings(),
                        RandomNumberEngine *random_number_engine = &random_global)
          : EnergyTermCommon(chain, "charmm36-cmap", settings, random_number_engine) {

          // Get CMAP data from the Gromacs code.
          this->cmapdata = charmm36_cmap::setup_cmap();
          this->cmap_pairs = topology::generate_cmap_pairs(this->chain);
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36Cmap(const TermCharmm36Cmap &other,
                        RandomNumberEngine *random_number_engine,
                        int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain),
            cmapdata(other.cmapdata),
            cmap_pairs(other.cmap_pairs) {}


     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double cmap_energy = 0.0;

          for (unsigned int i = 0; i < this->cmap_pairs.size(); i++) {

              const int residue_index = (this->cmap_pairs)[i].residue_index;
              const unsigned int cmap_type_index = (this->cmap_pairs)[i].cmap_type_index;
              // const unsigned int cmap_type_index = 0;

              const double phi = (*(this->chain))[residue_index].get_phi();
              const double psi = (*(this->chain))[residue_index].get_psi();

               cmap_energy += kernel::cmap_energy(phi, psi, cmap_type_index, this->cmapdata);

               // printf("             CMAP E = %12.4f kJ/mol\n", cmap_energy);
          }

          return cmap_energy / 4.184;
     }

};

}

#endif
