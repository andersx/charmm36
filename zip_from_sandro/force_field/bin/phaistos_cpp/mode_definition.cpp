
// Insert if not already present
if (std::find(phaistos_modes.begin(), phaistos_modes.end(), "mc-dynamics") == phaistos_modes.end())
     phaistos_modes.push_back("mc-dynamics");

// Module: Initialize mode definition (defined in includes.cpp)
module_force_field::ModeDefinitionInitialization(target, chain, atom_types_default, implicit_energies_allowed);
