function libDef = definelibdiv1d_extend(libDef,varargin)
%% DEFINELIBDIV1D About definelibdiv1d.mlx
% This file defines the MATLAB interface to the library |libdiv1d|.
% 
% Commented sections represent C++ functionality that MATLAB cannot automatically 
% define. To include functionality, uncomment a section and provide values for 
% <SHAPE>, <DIRECTION>, etc. For more information, see <matlab:helpview(fullfile(docroot,'matlab','helptargets.map'),'cpp_define_interface') 
% Define MATLAB Interface for C++ Library>.
%% Setup
% Do not edit this setup section.
% libDef = clibgen.LibraryDefinition("libdiv1dData.xml");
% %% OutputFolder and Libraries
% libDef.OutputFolder = "/home/unix/derks/Desktop/codes/div1d_library";
% libDef.Libraries = "/home/unix/derks/Desktop/codes/div1d_library/obj/libdiv1d.so";

% The Setup is generated based on the present paths and directories in div1d_library_interface.m

%% C++ function |initialize_div1d_settings_| with MATLAB name |clib.libdiv1d.initialize_div1d_settings_|
% C++ Signature: void initialize_div1d_settings_(double * floatinnum,int * intinnum,signed 
% char * loginnum,double * floatinphys,int * intinphys,int * call_from_extern)
initialize_div1d_settings_Definition = addFunction(libDef, ...
    "void initialize_div1d_settings_(double * floatinnum,int * intinnum,signed char * loginnum,double * floatinphys,int * intinphys,int * call_from_extern)", ...
    "MATLABName", "clib.libdiv1d.initialize_div1d_settings_", ...
    "Description", "clib.libdiv1d.initialize_div1d_settings_ Representation of C++ function initialize_div1d_settings_."); % Modify help description values as needed.
defineArgument(initialize_div1d_settings_Definition, "floatinnum", "double", "input", 90); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_settings_Definition, "intinnum", "int32", "input", 30); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
defineArgument(initialize_div1d_settings_Definition, "loginnum", "int8", "input", 10); % <MLTYPE> can be "clib.array.libdiv1d.SignedChar", or "int8"
defineArgument(initialize_div1d_settings_Definition, "floatinphys", "double", "input", 141); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_settings_Definition, "intinphys", "int32", "input", 31); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
defineArgument(initialize_div1d_settings_Definition, "call_from_extern", "int32", "input", 1); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
validate(initialize_div1d_settings_Definition);
%% C++ function |initialize_div1d_arrays_| with MATLAB name |clib.libdiv1d.initialize_div1d_arrays_|
% C++ Signature: void initialize_div1d_arrays_(double * E_density,double * E_velocity,double 
% * E_temperature,double * E_neutral,double * E_neutral_velocity,double * E_molecule,double 
% * E_x,double * E_xcb,double * E_delta_x,double * E_delta_xcb,double * E_B_field,double 
% * E_B_field_cb,double * E_B_trans,double * E_B_trans_cb,double * E_R_cc,double 
% * E_R_cb,double * E_Area_extern,double * E_sintheta_cc,double * E_sintheta_cb,double 
% * E_sol_width_pol,double * E_sol_width_pol_cb,double * E_volumes,double * E_gas_puff_profile,double 
% * E_core_source_profile_q,double * E_core_source_profile_n,int * init_grid_fortran,int 
% * init_prof_fortran,int * Nx,int * E_i_omp,int * E_i_Xpoint,int * E_i_baffle,double 
% * E_mid_point,double * E_X_omp)
initialize_div1d_arrays_Definition = addFunction(libDef, ...
    "void initialize_div1d_arrays_(double * E_density,double * E_velocity,double * E_temperature,double * E_neutral,double * E_neutral_velocity,double * E_molecule,double * E_x,double * E_xcb,double * E_delta_x,double * E_delta_xcb,double * E_B_field,double * E_B_field_cb,double * E_B_trans,double * E_B_trans_cb,double * E_R_cc,double * E_R_cb,double * E_Area_extern,double * E_sintheta_cc,double * E_sintheta_cb,double * E_sol_width_pol,double * E_sol_width_pol_cb,double * E_volumes,double * E_gas_puff_profile,double * E_core_source_profile_q,double * E_core_source_profile_n,int * init_grid_fortran,int * init_prof_fortran,int * Nx,int * E_i_omp,int * E_i_Xpoint,int * E_i_baffle,double * E_mid_point,double * E_X_omp)", ...
    "MATLABName", "clib.libdiv1d.initialize_div1d_arrays_", ...
    "Description", "clib.libdiv1d.initialize_div1d_arrays_ Representation of C++ function initialize_div1d_arrays_."); % Modify help description values as needed.
defineArgument(initialize_div1d_arrays_Definition, "E_density", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_velocity", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_temperature", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_neutral", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_neutral_velocity", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_molecule", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_x", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_xcb", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_delta_x", "double", "inputoutput", 199); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_delta_xcb", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_B_field", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_B_field_cb", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_B_trans", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_B_trans_cb", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_R_cc", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_R_cb", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_Area_extern", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_sintheta_cc", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_sintheta_cb", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_sol_width_pol", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_sol_width_pol_cb", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_volumes", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_gas_puff_profile", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_core_source_profile_q", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_core_source_profile_n", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "init_grid_fortran", "int32", "input", 1); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
defineArgument(initialize_div1d_arrays_Definition, "init_prof_fortran", "int32", "input", 1); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
defineArgument(initialize_div1d_arrays_Definition, "Nx", "int32", "input", 1); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
defineArgument(initialize_div1d_arrays_Definition, "E_i_omp", "int32", "inputoutput", 1); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
defineArgument(initialize_div1d_arrays_Definition, "E_i_Xpoint", "int32", "inputoutput", 2); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
defineArgument(initialize_div1d_arrays_Definition, "E_i_baffle", "int32", "inputoutput", 2); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
defineArgument(initialize_div1d_arrays_Definition, "E_mid_point", "double", "inputoutput",1); % <MLTYPE> can be "double", or "double"
defineArgument(initialize_div1d_arrays_Definition, "E_X_omp", "double", "inputoutput",1); % <MLTYPE> can be "double", or "double"
validate(initialize_div1d_arrays_Definition);
%% C++ function |run_div1d_| with MATLAB name |clib.libdiv1d.run_div1d_|
% C++ Signature: void run_div1d_(double * E_density,double * E_velocity,double 
% * E_temperature,double * E_neutral,double * E_neutral_velocity,double * E_molecule,double 
% * E_Gamma_n,double * E_Gamma_mom,double * E_q_parallel,double * E_Gamma_neutral,double 
% * E_Gamma_mom_neutral,double * E_Gamma_molecule,double * E_core2sol_flux,double 
% * E_core2sol_mol,double * E_sol2extern_ion_flux,double * E_Source_n,double * 
% E_Source_v,double * E_Source_Q,double * E_Source_neutral,double * E_Source_vn,double 
% * E_Source_molecule,double * E_extern_neutral_density,double * E_extern_neutral_flux,double 
% * E_sol2extern_flux,double * E_extern2sol_flux,double * E_tar2extern_flux,double 
% * E_extern2core_flux,double * E_sum_sol2extern_ion_flux,double * E_Source_extern,double 
% * E_neutral_pump,double * E_extern_molecule_density,double * E_extern_molecule_flux,double 
% * E_sol2extern_mol,double * E_extern2sol_mol,double * E_tar2extern_mol,double 
% * E_extern2core_mol,double * E_sum_sol2extern_ion_mol,double * E_Source_extern_mol,double 
% * E_molecule_pump,double * E_core_density,double * E_Gamma_core2sol,double * 
% E_sol2core_flux,double * E_sol2core_mol,double * E_Source_core,double * start_time,double 
% * end_time,int * nout,double * delta_t,double * E_imp_con,double * E_neu,size_t 
% len,double * E_dneu,double * E_nb,double * E_mb,double * E_gas,double * E_rec,double 
% * E_qpar_x,double * E_red_frc,double * E_Q_core,double * E_Gamma_core,double 
% * E_core_neutral_density,double * E_neutral_puff,double * E_molecule_puff,double 
% * E_core_fuelling)
run_div1d_Definition = addFunction(libDef, ...
    "void run_div1d_(double * E_density,double * E_velocity,double * E_temperature,double * E_neutral,double * E_neutral_velocity,double * E_molecule,double * E_Gamma_n,double * E_Gamma_mom,double * E_q_parallel,double * E_Gamma_neutral,double * E_Gamma_mom_neutral,double * E_Gamma_molecule,double * E_core2sol_flux,double * E_core2sol_mol,double * E_sol2extern_ion_flux,double * E_Source_n,double * E_Source_v,double * E_Source_Q,double * E_Source_neutral,double * E_Source_vn,double * E_Source_molecule,double * E_extern_neutral_density,double * E_extern_neutral_flux,double * E_sol2extern_flux,double * E_extern2sol_flux,double * E_tar2extern_flux,double * E_extern2core_flux,double * E_sum_sol2extern_ion_flux,double * E_Source_extern,double * E_neutral_pump,double * E_extern_molecule_density,double * E_extern_molecule_flux,double * E_sol2extern_mol,double * E_extern2sol_mol,double * E_tar2extern_mol,double * E_extern2core_mol,double * E_sum_sol2extern_ion_mol,double * E_Source_extern_mol,double * E_molecule_pump,double * E_core_density,double * E_Gamma_core2sol,double * E_sol2core_flux,double * E_sol2core_mol,double * E_Source_core,double * start_time,double * end_time,int * nout,double * delta_t,double * E_imp_con,double * E_neu,double * E_dneu,double * E_nb,double * E_mb,double * E_gas,double * E_rec,double * E_qpar_x,double * E_red_frc,double * E_Q_core,double * E_Gamma_core,double * E_core_neutral_density,double * E_neutral_puff,double * E_molecule_puff,double * E_core_fuelling)", ...
    "MATLABName", "clib.libdiv1d.run_div1d_", ...
    "Description", "clib.libdiv1d.run_div1d_ Representation of C++ function run_div1d_."); % Modify help description values as needed.
defineArgument(run_div1d_Definition, "E_density", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_velocity", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_temperature", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_neutral", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_neutral_velocity", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_molecule", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Gamma_n", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Gamma_mom", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_q_parallel", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Gamma_neutral", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Gamma_mom_neutral", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Gamma_molecule", "double", "inputoutput", 201); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_core2sol_flux", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_core2sol_mol", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_sol2extern_ion_flux", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Source_n", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Source_v", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Source_Q", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Source_neutral", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Source_vn", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Source_molecule", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_extern_neutral_density", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_extern_neutral_flux", "double", "inputoutput", 3); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_sol2extern_flux", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_extern2sol_flux", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_tar2extern_flux", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_extern2core_flux", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_sum_sol2extern_ion_flux", "double", "inputoutput", 1); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Source_extern", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_neutral_pump", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_extern_molecule_density", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_extern_molecule_flux", "double", "inputoutput", 3); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_sol2extern_mol", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_extern2sol_mol", "double", "inputoutput", 200); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_tar2extern_mol", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_extern2core_mol", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_sum_sol2extern_ion_mol", "double", "inputoutput", 1); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Source_extern_mol", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_molecule_pump", "double", "inputoutput", 5); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_core_density", "double", "inputoutput", 1); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Gamma_core2sol", "double", "inputoutput", 1); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_sol2core_flux", "double", "inputoutput", 1); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_sol2core_mol", "double", "inputoutput", 1); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Source_core", "double", "inputoutput", 1); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "start_time", "double", "input", 1); % <MLTYPE> TYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "end_time", "double", "inputoutput", 1); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "nout", "int32", "input",1); % <MLTYPE> can be "clib.array.libdiv1d.Int", or "int32"
defineArgument(run_div1d_Definition, "delta_t", "double", "input", 1); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_imp_con", "double", "input", [5,10]); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_neu", "double", "input", 10); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_dneu", "double", "input", 10); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_nb", "double", "input", [5,10]); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_mb", "double", "input", [5,10]); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_gas", "double", "input", 10); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_rec", "double", "input", 10); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_qpar_x", "double", "input", 10); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_red_frc", "double", "input", 10); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Q_core", "double", "input", 10); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_Gamma_core", "double", "input", 10); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_core_neutral_density", "double", "inputoutput", 10); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_neutral_puff", "double", "input",[5,10]); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_molecule_puff", "double", "input",[5,10]); % <MLTYPE> can be "double", or "double"
defineArgument(run_div1d_Definition, "E_core_fuelling", "double", "input",10); % <MLTYPE> can be "double", or "double"
validate(run_div1d_Definition);

%% Validate the library definition
validate(libDef);
end