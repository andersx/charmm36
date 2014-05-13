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
          option = options[prefix+"-gromacs-angle-bend"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermGromacsAngleBend::Settings Settings;

               // Add energy term
               energy->add_term(new TermGromacsAngleBend(chain,
                                                   options.get_settings<Settings>(option,i)));
          }

          // Bond Stretch
          option = options[prefix+"-gromacs-bond-stretch"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermGromacsBondStretch::Settings Settings;

               // Add energy term
               energy->add_term(new TermGromacsBondStretch(chain,
                                                   options.get_settings<Settings>(option,i)));
          }

          // Coulomb
          option = options[prefix+"-gromacs-coulomb"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermGromacsCoulomb::Settings Settings;

               // Add energy term
               energy->add_term(new TermGromacsCoulomb(chain,
                                                   options.get_settings<Settings>(option,i)));
          }

          // EEF1
          option = options[prefix+"-gromacs-eef1"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermGromacsEef1::Settings Settings;

               // Add energy term
               energy->add_term(new TermGromacsEef1(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Torsion
          option = options[prefix+"-gromacs-torsion"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermGromacsTorsion::Settings Settings;

               // Add energy term
               energy->add_term(new TermGromacsTorsion(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Improper torsion
          option = options[prefix+"-gromacs-imptor"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermGromacsImptor::Settings Settings;

               // Add energy term
               energy->add_term(new TermGromacsImptor(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
          // Lennard Jones
          option = options[prefix+"-gromacs-vdw"];
          for (int i=0; i<option.occurrences(); ++i) {

               // Settings typedef
               typedef TermGromacsVdw::Settings Settings;

               // Add energy term
               energy->add_term(new TermGromacsVdw(chain,
                                                   options.get_settings<Settings>(option,i)));
          }
     }

};

}
