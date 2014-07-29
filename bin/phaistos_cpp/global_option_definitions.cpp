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
          for (int counter = occurrences[prefix+"-charmm-angle-bend"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmmAngleBend EnergyTerm;
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
                         prefix+"-charmm-angle-bend", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Bond Stretch
          for (int counter = occurrences[prefix+"-charmm-bond-stretch"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmmBondStretch EnergyTerm;
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
                         prefix+"-charmm-bond-stretch", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // CMAP term
          for (int counter = occurrences[prefix+"-charmm-cmap"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmmCmap EnergyTerm;
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
                         prefix+"-charmm-cmap", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Coulomb
          for (int counter = occurrences[prefix+"-charmm-coulomb"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmmCoulomb EnergyTerm;
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
                         prefix+"-charmm-coulomb", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Implicit solvent
          for (int counter = occurrences[prefix+"-charmm-implicit-solvent"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmmImplicitSolvent EnergyTerm;
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
                         prefix+"-charmm-implicit-solvent", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Torsion
          for (int counter = occurrences[prefix+"-charmm-torsion"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmmTorsion EnergyTerm;
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
                         prefix+"-charmm-torsion", settings,
                         make_vector()),
                         super_group, counter==1);
          }

          // Improper torsion
          for (int counter = occurrences[prefix+"-charmm-improper-torsion"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmmImproperTorsion EnergyTerm;
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
                         prefix+"-charmm-improper-torsion", settings,
                         make_vector()),
                         super_group, counter==1);
          }

          // Lennard Jones
          for (int counter = occurrences[prefix+"-charmm-vdw"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmmVdw EnergyTerm;
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
                         prefix+"-charmm-vdw", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // // Non-Bonded
          for (int counter = occurrences[prefix+"-charmm-non-bonded"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmmNonBonded EnergyTerm;
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
                         "CHARMM36/EEF1-SB van der Waals, Coulomb and implicit solvent terms (" + prefix + ")",
                         prefix+"-charmm-non-bonded", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Bonded cached terms
          for (int counter = occurrences[prefix+"-charmm-bonded-cached"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmmBondedCached EnergyTerm;
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
                         prefix+"-charmm-bonded-cached", settings,
                         make_vector(
                             make_vector(std::string("ignore-bond-angles"),
                                         std::string("Ignore bond angle terms."),
                                          &settings->ignore_bond_angles),
                             make_vector(std::string("ignore-bond-stretch"),
                                         std::string("Ignore bond stretch terms."),
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
          for (int counter = occurrences[prefix+"-charmm-non-bonded-cached"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmmNonBondedCached EnergyTerm;
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
                         prefix+"-charmm-non-bonded-cached", settings,
                         make_vector()),
                    super_group, counter==1);
          }

     }


};

}
