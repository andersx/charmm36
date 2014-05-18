namespace module_force_field {

// Module: energy term initialization
struct EnergyInitialization {


     // Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyInitialization(const Options &options, CHAIN_TYPE *chain, DBN_TYPE *dbn,
                                       Energy<CHAIN_TYPE> *energy, std::vector<RandomNumberEngine *> *random_number_generators,
                                       std::string prefix="") {
     }

     // Constructor - template specific case
     template <typename DBN_TYPE>
     EnergyInitialization(const Options &options, ChainFB *chain, DBN_TYPE *dbn,
                                       Energy<ChainFB> *energy, std::vector<RandomNumberEngine *> *random_number_generators,
                                       std::string prefix="") {

          Options::OptionValue option;

          // Angle Bend
          option = options[prefix+"-charmm36-angle-bend"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmm36AngleBend::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmm36AngleBend(chain,
                                                   options.get_settings<Settings>(option,i)));
          }

          // Bond Stretch
          option = options[prefix+"-charmm36-bond-stretch"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmm36BondStretch::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmm36BondStretch(chain,
                                                   options.get_settings<Settings>(option,i)));
          }

          // Coulomb
          option = options[prefix+"-charmm36-coulomb"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmm36Coulomb::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmm36Coulomb(chain,
                                                   options.get_settings<Settings>(option,i)));
          }

          // EEF1
          option = options[prefix+"-charmm36-eef1"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmm36Eef1::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmm36Eef1(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Torsion
          option = options[prefix+"-charmm36-torsion"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmm36Torsion::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmm36Torsion(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Improper torsion
          option = options[prefix+"-charmm36-imptor"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmm36Imptor::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmm36Imptor(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Lennard Jones
          option = options[prefix+"-charmm36-vdw"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmm36Vdw::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmm36Vdw(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Non-bonded cached terms
          option = options[prefix+"-charmm36-non-bonded"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmm36NonBonded::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmm36NonBonded(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Non-bonded cached terms
          option = options[prefix+"-charmm36-non-bonded-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmm36NonBondedCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmm36NonBondedCached(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
     }

};

}
