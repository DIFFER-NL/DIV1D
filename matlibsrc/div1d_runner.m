function [out] = div1d_runner(in)
out = in;
% parses structs to the library and runs a timestep
[out.density, out.velocity, out.temperature, out.neutral_density, out.neutral_velocity, out.molecule, ... 6
 out.Gamma_n, out.Gamma_mom, out.q_parallel, out.Gamma_neutral, out.Gamma_mom_neutral, out.Gamma_molecule, ...12
 out.core2sol_flux, out.core2sol_mol, out.sol2extern_ion_flux,  ...
 out.Source_n, out.Source_v, out.Source_Q, out.Source_neutral, out.Source_vn, out.Source_molecule,  ...21
 out.extern_neutral_density,  ... 22
 out.extern_neutral_flux, out.sol2extern_flux, out.extern2sol_flux, out.tar2extern_flux,  ... 26
 out.extern2core_flux, out.sum_sol2extern_ion_flux, ... 28
 out.Source_extern, out.neutral_pump, ...30
 out.extern_molecule_density,  ...31
 out.extern_molecule_flux, out.sol2extern_mol, out.extern2sol_mol, out.tar2extern_mol, ...35
 out.extern2core_mol, out.sum_sol2extern_ion_mol,  ...37
 out.Source_extern_mol, out.molecule_pump, ...39
 out.core_density, out.Gamma_core2sol, out.sol2core_flux, out.sol2core_mol, out.Source_core, ...
 out.end_time, out.core_neutral_density]...
 = clib.libdiv1d.run_div1d_(	...
 in.density, in.velocity, in.temperature, in.neutral_density, in.neutral_velocity, in.molecule, ... 6
 in.Gamma_n, in.Gamma_mom, in.q_parallel, in.Gamma_neutral, in.Gamma_mom_neutral, in.Gamma_molecule, ...12
 in.core2sol_flux, in.core2sol_mol, in.sol2extern_ion_flux,  ...
 in.Source_n, in.Source_v, in.Source_Q, in.Source_neutral, in.Source_vn, in.Source_molecule,  ...21
 in.extern_neutral_density,  ... 22
 in.extern_neutral_flux, in.sol2extern_flux, in.extern2sol_flux, in.tar2extern_flux,  ... 26
 in.extern2core_flux, in.sum_sol2extern_ion_flux, ... 28
 in.Source_extern, in.neutral_pump, ...30
 in.extern_molecule_density,  ...31
 in.extern_molecule_flux, in.sol2extern_mol, in.extern2sol_mol, in.tar2extern_mol, ...35
 in.extern2core_mol, in.sum_sol2extern_ion_mol,  ...37
 in.Source_extern_mol, in.molecule_pump, ...39
 in.core_density, in.Gamma_core2sol, in.sol2core_flux, in.sol2core_mol, in.Source_core, ...44
 in.start_time, in.end_time, in.nout, in.delta_t,  ...48
 in.imp_con, in.neu, in.dneu, in.nb, in.mb, ... 53 
 in.gas, in.rec, in.qpar_x, in.red_frc, ... 57
 in.Q_core, in.Gamma_core, in.core_neutral_density , ...60
 in.neutral_puff, in.molecule_puff, in.core_fuelling); ...63

out.time = out.start_time:in.delta_t*in.nout:out.end_time-in.delta_t*in.nout;
end