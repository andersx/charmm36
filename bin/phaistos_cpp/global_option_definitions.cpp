namespace module_force_field {

//! Module: energy term initialization
template <typename SETTINGS_MODIFIER>
struct EnergyOptions {

     //! Constructor - general case: do nothing
     template <typename CHAIN_TYPE, typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   CHAIN_TYPE *chain,
                   DBN_TYPE *dbn) {
     }

     //! Constructor - ChainFB specific case
     template <typename DBN_TYPE>
     EnergyOptions(ProgramOptionParser &target,
                   const ProgramOptionParser::Filter &occurrences,
                   std::string super_group,
                   std::string prefix,
                   ChainFB *chain,
                   DBN_TYPE *dbn) {

          // Facilitate use of make_vector
          using namespace boost::fusion;

          // Define defaults for the different modes
          ModeDefinitions mode_definitions(target, chain);

          // Angle Bend
          for (int counter = occurrences[prefix+"-charmm36-angle-bend"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36AngleBend EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "CHARMM36/EEF1-SB angle bend term (" + prefix + ")",
                         prefix+"-charmm36-angle-bend", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Bond Stretch
          for (int counter = occurrences[prefix+"-charmm36-bond-stretch"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36BondStretch EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "CHARMM36/EEF1-SB bond stretch term (" + prefix + ")",
                         prefix+"-charmm36-bond-stretch", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // CMAP term
          for (int counter = occurrences[prefix+"-charmm36-cmap"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36Cmap EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "CHARMM36/EEF1-SB CMAP correction term (" + prefix + ")",
                         prefix+"-charmm36-cmap", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Coulomb
          for (int counter = occurrences[prefix+"-charmm36-coulomb"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36Coulomb EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "CHARMM36/EEF1-SB Coulomb term (" + prefix + ")",
                         prefix+"-charmm36-coulomb", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Implicit solvent
          for (int counter = occurrences[prefix+"-charmm36-implicit-solvent"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36ImplicitSolvent EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "CHARMM36/EEF1-SB implicit solvation term (" + prefix + ")",
                         prefix+"-charmm36-implicit-solvent", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Torsion
          for (int counter = occurrences[prefix+"-charmm36-torsion"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36Torsion EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "CHARMM36/EEF1-SB torsion angle term (" + prefix + ")",
                         prefix+"-charmm36-torsion", settings,
                         make_vector()),
                         super_group, counter==1);
          }

          // Improper torsion
          for (int counter = occurrences[prefix+"-charmm36-improper-torsion"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36ImproperTorsion EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "CHARMM36/EEF1-SB improper torsion angle term (" + prefix + ")",
                         prefix+"-charmm36-imptroper", settings,
                         make_vector()),
                         super_group, counter==1);
          }

          // Lennard Jones
          for (int counter = occurrences[prefix+"-charmm36-vdw"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36Vdw EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "CHARMM36/EEF1-SB van der Waals term (" + prefix + ")",
                         prefix+"-charmm36-vdw", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // // Non-Bonded
          // for (int counter = occurrences[prefix+"-charmm36-non-bonded"]; counter > 0; counter--) {

          //      // Create settings object
          //      typedef TermCharmm36NonBonded EnergyTerm;
          //      typedef EnergyTerm::Settings Settings;
          //      boost::shared_ptr<Settings> settings(
          //           SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

          //      // If temperature has been specified, set weight based on that value
          //      if(target.has_key("temperature")) {
          //           settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
          //      }

          //      // Add options
          //      target.add(
          //           target.create_options(
          //                DefineEnergyCommonOptions(),
          //                "Van der Waal + Coulomb + EEF1-SB terms (" + prefix + ")",
          //                prefix+"-charmm36-non-bonded", settings,
          //                make_vector()),
          //           super_group, counter==1);
          // }

          // Bonded cached terms
          for (int counter = occurrences[prefix+"-charmm36-bonded-cached"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36BondedCached EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "Cached CHARMM36/EEF1-SB angle bend, bond stretch, torsion angle, improper torsion angles and CMAP correction terms (" + prefix + ")",
                         prefix+"-charmm36-bonded-cached", settings,
                         make_vector(
                             make_vector(std::string("ignore-bond-angles"),
                                         std::string("Ignore bond angle terms."),
                                          &settings->ignore_bond_angles),
                             make_vector(std::string("ignore-bond-stretch"),
                                         std::string("Ignore bond stretch terms (Note: Ignored by default)."),
                                          &settings->ignore_bond_stretch),
                             make_vector(std::string("ignore-torsion-angles"),
                                         std::string("Ignore torsion angle terms."),
                                          &settings->ignore_torsion_angles),
                             make_vector(std::string("ignore-improper-torsion-angles"),
                                         std::string("Ignore improper torsion angle terms."),
                                          &settings->ignore_improper_torsion_angles),
                             make_vector(std::string("ignore-cmap-correction"),
                                         std::string("Ignore CMAP correction terms."),
                                          &settings->ignore_cmap_correction)
                        )),
                    super_group, counter==1);
          }

          // Non-Bonded chaced
          for (int counter = occurrences[prefix+"-charmm36-non-bonded-cached"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36NonBondedCached EnergyTerm;
               typedef EnergyTerm::Settings Settings;
               boost::shared_ptr<Settings> settings(
                    SETTINGS_MODIFIER().template modify<EnergyTerm>(new Settings(), prefix));

               // If temperature has been specified, set weight based on that value
               if(target.has_key("temperature")) {
                    settings->weight = temperature_to_one_over_k(target["temperature"].as<double>());
               }

               // Add options
               target.add(
                    target.create_options(
                         DefineEnergyCommonOptions(),
                         "Cached CHARMM36/EEF1-SB van der Waals, Coulomb and implicit solvent terms (" + prefix + ")",
                         prefix+"-charmm36-non-bonded-cached", settings,
                         make_vector()),
                    super_group, counter==1);
          }

     }


};

}
