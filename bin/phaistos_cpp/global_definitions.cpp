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
          option = options[prefix+"-charmm-angle-bend"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmmAngleBend::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmmAngleBend(chain,
                                                   options.get_settings<Settings>(option,i)));
          }

          // Bond Stretch
          option = options[prefix+"-charmm-bond-stretch"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmmBondStretch::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmmBondStretch(chain,
                                                   options.get_settings<Settings>(option,i)));
          }

          // Coulomb
          option = options[prefix+"-charmm-coulomb"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmmCoulomb::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmmCoulomb(chain,
                                                   options.get_settings<Settings>(option,i)));
          }

          // Implicit solvent
          option = options[prefix+"-charmm-implicit-solvent"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmmImplicitSolvent::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmmImplicitSolvent(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Torsion
          option = options[prefix+"-charmm-torsion"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmmTorsion::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmmTorsion(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Improper torsion
          option = options[prefix+"-charmm-ImproperTorsion"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmmImproperTorsion::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmmImproperTorsion(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Lennard Jones
          option = options[prefix+"-charmm-vdw"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmmVdw::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmmVdw(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // // Non-bonded terms
          option = options[prefix+"-charmm-non-bonded"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmmNonBonded::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmmNonBonded(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Non-bonded cached terms
          option = options[prefix+"-charmm-non-bonded-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmmNonBondedCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmmNonBondedCached(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // CMAP term
          option = options[prefix+"-charmm-cmap"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmmCmap::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmmCmap(chain,
                                                options.get_settings<Settings>(option,i)));
          }
          // bonded-cached term
          option = options[prefix+"-charmm-bonded-cached"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermCharmmBondedCached::Settings Settings;

               // Add energy term
               energy->add_term(new TermCharmmBondedCached(chain,
                                                options.get_settings<Settings>(option,i)));
          }
     }

};

}
