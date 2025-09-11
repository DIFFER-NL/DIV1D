%% About definelibdiv1d.mlx
% This file defines the MATLAB interface to the library |libdiv1d|.
%
% Commented sections represent C++ functionality that MATLAB cannot automatically define. To include
% functionality, uncomment a section and provide values for &lt;SHAPE&gt;, &lt;DIRECTION&gt;, etc. For more
% information, see <matlab:helpview(fullfile(docroot,'matlab','helptargets.map'),'cpp_define_interface') Define MATLAB Interface for C++ Library>.



%% Setup
% Do not edit this setup section.
function libDef = definelibdiv1d()
libDef = clibgen.LibraryDefinition("libdiv1dData.xml");
%% OutputFolder and Libraries 
libDef.OutputFolder = "/home/unix/derks/Desktop/projects/dynamics/models/div1d";
libDef.Libraries = "/home/unix/derks/Desktop/projects/dynamics/models/div1d/obj/libdiv1d.so";

%% C++ function |initialize_div1d_settings_| with MATLAB name |clib.libdiv1d.initialize_div1d_settings_|
% C++ Signature: void initialize_div1d_settings_(double * floatinnum,int * intinnum,signed char * loginnum,double * floatinphys,int * intinphys,int * call_from_extern)
%initialize_div1d_settings_Definition = addFunction(libDef, ...
%    "void initialize_div1d_settings_(double * floatinnum,int * intinnum,signed char * loginnum,double * floatinphys,int * intinphys,int * call_from_extern)", ...
%    "MATLABName", "clib.libdiv1d.initialize_div1d_settings_", ...
%    "Description", "clib.libdiv1d.initialize_div1d_settings_ Representation of C++ function initialize_div1d_settings_."); % Modify help description values as needed.
%defineArgument(initialize_div1d_settings_Definition, "floatinnum", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_settings_Definition, "intinnum", "clib.array.libdiv1d.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
%defineArgument(initialize_div1d_settings_Definition, "loginnum", "clib.array.libdiv1d.SignedChar", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.SignedChar", or "int8"
%defineArgument(initialize_div1d_settings_Definition, "floatinphys", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_settings_Definition, "intinphys", "clib.array.libdiv1d.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
%defineArgument(initialize_div1d_settings_Definition, "call_from_extern", "clib.array.libdiv1d.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
%validate(initialize_div1d_settings_Definition);

%% C++ function |initialize_div1d_arrays_| with MATLAB name |clib.libdiv1d.initialize_div1d_arrays_|
% C++ Signature: void initialize_div1d_arrays_(double * E_density,double * E_velocity,double * E_temperature,double * E_neutral,double * E_neutral_velocity,double * E_molecule,double * E_x,double * E_xcb,double * E_delta_x,double * E_delta_xcb,double * E_B_field,double * E_B_field_cb,double * E_B_trans,double * E_B_trans_cb,double * E_R_cc,double * E_R_cb,double * E_Area_extern,double * E_sintheta_cc,double * E_sintheta_cb,double * E_sol_width_pol,double * E_sol_width_pol_cb,double * E_volumes,double * E_gas_puff_profile,double * E_core_source_profile_q,double * E_core_source_profile_n,int * init_grid_fortran,int * init_prof_fortran,int * Nx,int * E_i_omp,int * E_i_Xpoint,int * E_i_baffle,double * E_mid_point,double * E_X_omp)
%initialize_div1d_arrays_Definition = addFunction(libDef, ...
%    "void initialize_div1d_arrays_(double * E_density,double * E_velocity,double * E_temperature,double * E_neutral,double * E_neutral_velocity,double * E_molecule,double * E_x,double * E_xcb,double * E_delta_x,double * E_delta_xcb,double * E_B_field,double * E_B_field_cb,double * E_B_trans,double * E_B_trans_cb,double * E_R_cc,double * E_R_cb,double * E_Area_extern,double * E_sintheta_cc,double * E_sintheta_cb,double * E_sol_width_pol,double * E_sol_width_pol_cb,double * E_volumes,double * E_gas_puff_profile,double * E_core_source_profile_q,double * E_core_source_profile_n,int * init_grid_fortran,int * init_prof_fortran,int * Nx,int * E_i_omp,int * E_i_Xpoint,int * E_i_baffle,double * E_mid_point,double * E_X_omp)", ...
%    "MATLABName", "clib.libdiv1d.initialize_div1d_arrays_", ...
%    "Description", "clib.libdiv1d.initialize_div1d_arrays_ Representation of C++ function initialize_div1d_arrays_."); % Modify help description values as needed.
%defineArgument(initialize_div1d_arrays_Definition, "E_density", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_velocity", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_temperature", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_neutral", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_neutral_velocity", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_molecule", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_x", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_xcb", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_delta_x", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_delta_xcb", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_B_field", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_B_field_cb", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_B_trans", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_B_trans_cb", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_R_cc", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_R_cb", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_Area_extern", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_sintheta_cc", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_sintheta_cb", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_sol_width_pol", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_sol_width_pol_cb", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_volumes", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_gas_puff_profile", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_core_source_profile_q", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_core_source_profile_n", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "init_grid_fortran", "clib.array.libdiv1d.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
%defineArgument(initialize_div1d_arrays_Definition, "init_prof_fortran", "clib.array.libdiv1d.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
%defineArgument(initialize_div1d_arrays_Definition, "Nx", "clib.array.libdiv1d.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
%defineArgument(initialize_div1d_arrays_Definition, "E_i_omp", "clib.array.libdiv1d.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
%defineArgument(initialize_div1d_arrays_Definition, "E_i_Xpoint", "clib.array.libdiv1d.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
%defineArgument(initialize_div1d_arrays_Definition, "E_i_baffle", "clib.array.libdiv1d.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
%defineArgument(initialize_div1d_arrays_Definition, "E_mid_point", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(initialize_div1d_arrays_Definition, "E_X_omp", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%validate(initialize_div1d_arrays_Definition);

%% C++ function |run_div1d_| with MATLAB name |clib.libdiv1d.run_div1d_|
% C++ Signature: void run_div1d_(double * E_density,double * E_velocity,double * E_temperature,double * E_neutral,double * E_neutral_velocity,double * E_molecule,double * E_Gamma_n,double * E_Gamma_mom,double * E_q_parallel,double * E_Gamma_neutral,double * E_Gamma_mom_neutral,double * E_Gamma_molecule,double * E_core2sol_flux,double * E_core2sol_mol,double * E_sol2extern_ion_flux,double * E_Source_n,double * E_Source_v,double * E_Source_Q,double * E_Source_neutral,double * E_Source_vn,double * E_Source_molecule,double * E_extern_neutral_density,double * E_extern_neutral_flux,double * E_sol2extern_flux,double * E_extern2sol_flux,double * E_tar2extern_flux,double * E_extern2core_flux,double * E_sum_sol2extern_ion_flux,double * E_Source_extern,double * E_neutral_pump,double * E_extern_molecule_density,double * E_extern_molecule_flux,double * E_sol2extern_mol,double * E_extern2sol_mol,double * E_tar2extern_mol,double * E_extern2core_mol,double * E_sum_sol2extern_ion_mol,double * E_Source_extern_mol,double * E_molecule_pump,double * E_core_density,double * E_Gamma_core2sol,double * E_sol2core_flux,double * E_sol2core_mol,double * E_Source_core,double * start_time,double * end_time,int * nout,double * delta_t,double * E_imp_con,double * E_neu,double * E_dneu,double * E_nb,double * E_mb,double * E_gas,double * E_rec,double * E_qpar_x,double * E_red_frc,double * E_Q_core,double * E_Gamma_core,double * E_core_neutral_density,double * E_neutral_puff,double * E_molecule_puff,double * E_core_fuelling)
%run_div1d_Definition = addFunction(libDef, ...
%    "void run_div1d_(double * E_density,double * E_velocity,double * E_temperature,double * E_neutral,double * E_neutral_velocity,double * E_molecule,double * E_Gamma_n,double * E_Gamma_mom,double * E_q_parallel,double * E_Gamma_neutral,double * E_Gamma_mom_neutral,double * E_Gamma_molecule,double * E_core2sol_flux,double * E_core2sol_mol,double * E_sol2extern_ion_flux,double * E_Source_n,double * E_Source_v,double * E_Source_Q,double * E_Source_neutral,double * E_Source_vn,double * E_Source_molecule,double * E_extern_neutral_density,double * E_extern_neutral_flux,double * E_sol2extern_flux,double * E_extern2sol_flux,double * E_tar2extern_flux,double * E_extern2core_flux,double * E_sum_sol2extern_ion_flux,double * E_Source_extern,double * E_neutral_pump,double * E_extern_molecule_density,double * E_extern_molecule_flux,double * E_sol2extern_mol,double * E_extern2sol_mol,double * E_tar2extern_mol,double * E_extern2core_mol,double * E_sum_sol2extern_ion_mol,double * E_Source_extern_mol,double * E_molecule_pump,double * E_core_density,double * E_Gamma_core2sol,double * E_sol2core_flux,double * E_sol2core_mol,double * E_Source_core,double * start_time,double * end_time,int * nout,double * delta_t,double * E_imp_con,double * E_neu,double * E_dneu,double * E_nb,double * E_mb,double * E_gas,double * E_rec,double * E_qpar_x,double * E_red_frc,double * E_Q_core,double * E_Gamma_core,double * E_core_neutral_density,double * E_neutral_puff,double * E_molecule_puff,double * E_core_fuelling)", ...
%    "MATLABName", "clib.libdiv1d.run_div1d_", ...
%    "Description", "clib.libdiv1d.run_div1d_ Representation of C++ function run_div1d_."); % Modify help description values as needed.
%defineArgument(run_div1d_Definition, "E_density", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_velocity", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_temperature", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_neutral", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_neutral_velocity", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_molecule", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Gamma_n", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Gamma_mom", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_q_parallel", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Gamma_neutral", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Gamma_mom_neutral", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Gamma_molecule", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_core2sol_flux", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_core2sol_mol", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_sol2extern_ion_flux", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Source_n", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Source_v", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Source_Q", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Source_neutral", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Source_vn", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Source_molecule", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_extern_neutral_density", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_extern_neutral_flux", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_sol2extern_flux", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_extern2sol_flux", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_tar2extern_flux", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_extern2core_flux", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_sum_sol2extern_ion_flux", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Source_extern", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_neutral_pump", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_extern_molecule_density", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_extern_molecule_flux", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_sol2extern_mol", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_extern2sol_mol", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_tar2extern_mol", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_extern2core_mol", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_sum_sol2extern_ion_mol", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Source_extern_mol", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_molecule_pump", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_core_density", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Gamma_core2sol", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_sol2core_flux", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_sol2core_mol", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Source_core", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "start_time", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "end_time", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "nout", "clib.array.libdiv1d.Int", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
%defineArgument(run_div1d_Definition, "delta_t", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_imp_con", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_neu", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_dneu", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_nb", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_mb", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_gas", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_rec", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_qpar_x", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_red_frc", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Q_core", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_Gamma_core", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_core_neutral_density", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_neutral_puff", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_molecule_puff", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%defineArgument(run_div1d_Definition, "E_core_fuelling", "clib.array.libdiv1d.Double", "input", <SHAPE>); % <MLTYPE> can be "clib.array.libdiv1d.Double", or "double"
%validate(run_div1d_Definition);

%% Validate the library definition
validate(libDef);

end
