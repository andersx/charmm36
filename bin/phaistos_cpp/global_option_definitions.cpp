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
          for (int counter = occurrences[prefix+"-gromacs-angle-bend"]; counter > 0; counter--) {

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
                         "Gromacs angle bend term (" + prefix + ")",
                         prefix+"-gromacs-angle-bend", settings,
                         make_vector(
                             make_vector(std::string("omit-sidechains"),
                             std::string("Whether to omit sidechains??"),
                             &settings->omit_sidechains)

                        )),
                    super_group, counter==1);
          }

          // Bond Stretch
          for (int counter = occurrences[prefix+"-gromacs-bond-stretch"]; counter > 0; counter--) {

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
                         "Gromacs bond stretch term (" + prefix + ")",
                         prefix+"-gromacs-bond-stretch", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Coulomb
          for (int counter = occurrences[prefix+"-gromacs-coulomb"]; counter > 0; counter--) {

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
                         "Gromacs Coulomb term (" + prefix + ")",
                         prefix+"-gromacs-coulomb", settings,
                         make_vector(
                             make_vector(std::string("e14factor"),
                                         std::string("E14FAC-coulomb"),
                                          &settings->E14FAC),
                             make_vector(std::string("epsilon"),
                                         std::string("Dielectric constant"),
                                         &settings->EPS),
                             make_vector(std::string("rdie"),
                                         std::string("Distance-dependent-dielectric??"),
                                         &settings->RDIE)
                        )),
                    super_group, counter==1);
          }

          // Coulomb
          for (int counter = occurrences[prefix+"-gromacs-eef1"]; counter > 0; counter--) {

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
                         "Gromacs EEF1 solvation term (" + prefix + ")",
                         prefix+"-gromacs-eef1", settings,
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
          for (int counter = occurrences[prefix+"-gromacs-torsion"]; counter > 0; counter--) {

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
                         "Gromacs torsion term (" + prefix + ")",
                         prefix+"-gromacs-torsion", settings,
                         make_vector()),
                         super_group, counter==1);
          }

          // Improper torsion
          for (int counter = occurrences[prefix+"-gromacs-imptor"]; counter > 0; counter--) {

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
                         "Gromacs improper torsion term (" + prefix + ")",
                         prefix+"-gromacs-imptor", settings,
                         make_vector()),
                         super_group, counter==1);
          }

          // Lennard Jones
          for (int counter = occurrences[prefix+"-gromacs-vdw"]; counter > 0; counter--) {

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
                         "Gromacs van der Waals term (" + prefix + ")",
                         prefix+"-gromacs-vdw", settings,
                         make_vector()),
                    super_group, counter==1);
          }

          // Non-Bonded chaced
          for (int counter = occurrences[prefix+"-gromacs-non-bonded-cached"]; counter > 0; counter--) {

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
                         prefix+"-gromacs-non-bonded-cached", settings,
                         make_vector()),
                    super_group, counter==1);
          }


     }


};

}
