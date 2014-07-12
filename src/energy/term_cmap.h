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
#include "term_cmap_tables.h"

namespace phaistos {

//! CMAP energy term
class TermCharmm36Cmap: public EnergyTermCommon<TermCharmm36Cmap, ChainFB> {

protected:

     //! For convenience, define local EnergyTermCommon
     typedef phaistos::EnergyTermCommon<TermCharmm36Cmap, ChainFB> EnergyTermCommon;

public:


     //! Table which contains the CMAP correction tables
     //! in Gromacs' internal formatting
     std::vector<std::vector<double> > cmap_data;

     //! Vector containing all terms in the CMAP correction.
     std::vector<topology::CmapInteraction> cmap_interactions;

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
          this->cmap_data = charmm36_cmap::setup_cmap();

          // Make a list of all CMAP interactions
          this->cmap_interactions = topology::generate_cmap_interactions(this->chain);
     }

     //! Copy constructor.
     //! \param other Source object from which copy is made
     //! \param random_number_engine Object from which random number generators can be created.
     //! \param thread_index Index indicating in which thread|rank the copy exists
     //! \param chain Molecule chain
     TermCharmm36Cmap(const TermCharmm36Cmap &other,
                        RandomNumberEngine *random_number_engine,
                        int thread_index, ChainFB *chain)
          : EnergyTermCommon(other, random_number_engine, thread_index, chain) {

          // Get CMAP data from the Gromacs code.
          this->cmap_data = charmm36_cmap::setup_cmap();

          // Make a list of all CMAP interactions
          this->cmap_interactions = topology::generate_cmap_interactions(this->chain);

     }

     //! Evaluate chain energy
     //! \param move_info object containing information about last move
     //! \return torsional potential energy of the chain in the object
     double evaluate(MoveInfo *move_info=NULL) {

          double cmap_energy = 0.0;

          for (unsigned int i = 0; i < this->cmap_interactions.size(); i++) {

               const int residue_index = (this->cmap_interactions)[i].residue_index;
               const unsigned int cmap_type_index = (this->cmap_interactions)[i].cmap_type_index;
               const double phi = (*(this->chain))[residue_index].get_phi();
               const double psi = (*(this->chain))[residue_index].get_psi();

               cmap_energy += charmm36_cmap::cmap_energy(phi, psi, cmap_type_index, this->cmap_data);

          }
          
          if (settings.debug > 0) {
               printf("             CMAP E = %15.6f kJ/mol\n", cmap_energy);
               printf("             CMAP E = %15.6f kcal/mol\n", cmap_energy * charmm36_constants::KJ_TO_KCAL);
          }

          return cmap_energy  * charmm36_constants::KJ_TO_KCAL;
     }

};

}

#endif
