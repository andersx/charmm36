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
                         "Charmm36 angle bend term (" + prefix + ")",
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
                         "Charmm36 bond stretch term (" + prefix + ")",
                         prefix+"-charmm36-bond-stretch", settings,
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
                         "Charmm36 Coulomb term (" + prefix + ")",
                         prefix+"-charmm36-coulomb", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Coulomb
          for (int counter = occurrences[prefix+"-charmm36-eef1"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36Eef1 EnergyTerm;
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
                         "Charmm36 EEF1 solvation term (" + prefix + ")",
                         prefix+"-charmm36-eef1", settings,
                         make_vector(
                             make_vector(std::string("solvation-filename"),
                                         std::string("solvation-filename"),
                                          &settings->solvation_filename),
                             make_vector(std::string("temperature"),
                                         std::string("Simulation Temperature."),
                                         &settings->temp)
                        )),
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
                         "Charmm36 torsion term (" + prefix + ")",
                         prefix+"-charmm36-torsion", settings,
                         make_vector()),
                         super_group, counter==1);
          }

          // Improper torsion
          for (int counter = occurrences[prefix+"-charmm36-imptor"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36Imptor EnergyTerm;
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
                         "Charmm36 improper torsion term (" + prefix + ")",
                         prefix+"-charmm36-imptor", settings,
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
                         "Charmm36 van der Waals term (" + prefix + ")",
                         prefix+"-charmm36-vdw", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Non-Bonded
          for (int counter = occurrences[prefix+"-charmm36-non-bonded"]; counter > 0; counter--) {

               // Create settings object
               typedef TermCharmm36NonBonded EnergyTerm;
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
                         "Van der Waal + Coulomb + EEF1-SB terms (" + prefix + ")",
                         prefix+"-charmm36-non-bonded", settings,
                         make_vector()),
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
                         "Cached van der Waal + Coulomb + EEF1-SB terms (" + prefix + ")",
                         prefix+"-charmm36-non-bonded-cached", settings,
                         make_vector()),
                    super_group, counter==1);
          }


     }


};

}
