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
          
          // Lennard Jones
          for (int counter = occurrences[prefix+"-vdw"]; counter > 0; counter--) {
          
               // Create settings object
               typedef TermVdw EnergyTerm;
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
                         "vdw: van der Waals term (" + prefix + ")",
                         prefix+"-vdw", settings,
                         make_vector()),
                    super_group, counter==1);
          }
          

     }

          
};

}
