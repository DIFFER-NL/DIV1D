% build the C interface to the DIV1D library

% run this script from the div1d folder
clear vars;
system('./build_libdiv1d.sh');
filePath = pwd;
addpath('obj');
addpath('src');
addpath('matlibsrc')
addpath('namelist')

hFile = strcat(filePath,'/src/libdiv1d.h');
libFile = strcat(filePath,'/obj/libdiv1d.so');
 try
clibgen.generateLibraryDefinition(hFile,Libraries=libFile,...
    CLinkage=true,verbose=true,OverwriteExistingDefinitionFiles=true) 

catch
disp('M definelibdivd.m(lx) already defined, moving on');
 end
% now go into .m and .mlx file to inform about shapes and sizes
%%
libpath = '/home/unix/derks/Desktop/projects/dynamics/models/div1d';
libDef = definelibdiv1d();
% Rsummary(definelibdiv1d(libpath))
summary(libDef)
libDef = definelibdiv1d_extend(libDef);
summary(libDef)

%% add libdiv1d.so to runtime path
% libpath = '/home/unix/derks/Desktop/projects/dynamics/models/div1d-debug/div1d';
% addpath('/home/unix/derks/Desktop/codes/div1d_library/obj')
div1dlibpath = [libpath,'/obj'];
addpath(div1dlibpath)

dllPath = 'rtPath'; 
%div1dlibpath = '/home/unix/derks/Desktop/codes/div1d_library/obj';
syspath = getenv('PATH'); 
setenv('PATH',[dllPath pathsep syspath div1dlibpath]);

%% now build it.
init_div1d = 0;
try 
% build(definelibdiv1d(libpath))
build(libDef)
addpath([libpath,'/libdiv1d'])
system('cp ./obj/libdiv1d.so ./libdiv1d/');
 % probably you have also initialized it

 init_div1d_settings = 0;
 init_div1d_arrays = 0;
catch
    init_div1d =0;
error( 'library already build or libdiv1d.so not finable')
end

% see if all libraries can be found
addpath( libpath,'/libdiv1d')
system('ldd libdiv1d/libdiv1dInterface.so');
%% define DIV1D object 
% settings
%addpath('/home/unix/derks/Desktop/codes/div1d_library/libdiv1d')

i = div1d_lib_inputs;
i.call_from_extern = 1;
%
if init_div1d_settings ==0
clib.libdiv1d.initialize_div1d_settings_( ...
            i.floatinnum, i.intinnum, i.loginnum, ...
            i.floatinphys,i.intinphys, i.call_from_extern  );
init_div1d_settings = 1
 else
 disp('M cannot initialize div1d twice, segmentation fault for second allocations')
 end

% works! 
%%
if init_div1d_arrays == 0
a = div1d_lib_arrays(i.intinnum(1)); 
[a.density, a.velocity, a.temperature, a.neutral_density, a.neutral_velocity, a.molecule, ...
		           a.x, a.xcb, a.delta_x, a.delta_xcb, a.b_fieldd, a.b_field_cb, a.b_trans, a.b_trans_cb, ...
			   a.r_cc, a.r_cb, a.area_extern, a.sintheta_cc, a.sintheta_cb, a.sol_width_pol, a.sol_width_pol_cb, a.volumes,...
			   a.gas_puff_profile, a.core_source_profile_q, a.core_source_profile_n,...
		       a.i_omp, a.i_xpoint, a.i_baffle, a.mid_point, a.x_omp] ...
    =clib.libdiv1d.initialize_div1d_arrays_(a.density, a.velocity, a.temperature, a.neutral_density, a.neutral_velocity, a.molecule, ...
		           a.x, a.xcb, a.delta_x, a.delta_xcb, a.b_field, a.b_field_cb, a.b_trans, a.b_trans_cb, ...
			   a.r_cc, a.r_cb, a.area_extern, a.sintheta_cc, a.sintheta_cb, a.sol_width_pol, a.sol_width_pol_cb, a.volumes,...
			   a.gas_puff_profile, a.core_source_profile_q, a.core_source_profile_n,...
		           a.init_grid_fortran, a.init_prof_fortran, i.intinnum(1), a.i_omp, a.i_xpoint, a.i_baffle, a.mid_point, a.x_omp);
 init_div1d_arrays = 1;
 else
 disp('M cannot setup div1d arrays twice, segmentation fault for second allocations')
 end
% works! but be carefull with i.intinnum(1)=Nx being the same here as in div1d_settings
%for
%while
%inputs(t) = 
%run_div1d
%run_jan
%end

%% now all has been initialized we can get the div1d input struct from the
% echo file.
echo_path = '/home/unix/derks/Desktop/projects/dynamics/models/div1d-debug/div1d/echo_div1d_inputs.txt';
input_struct = read_echo_div1d_inputs(echo_path);
input_struct.grid = a;

%%
r = div1d_runstruct(a.nx,i.intinnum(2), i);
r.density = a.density;
r.velocity = a.velocity;
r.neutral_density = a.neutral_density;
r.neutral_velocity = a.neutral_velocity;
r.molecule = a.molecule;


[r.density, r.velocity, r.temperature, r.neutral_density, r.neutral_velocity, r.molecule, ... 6
 r.Gamma_n, r.Gamma_mom, r.q_parallel, r.Gamma_neutral, r.Gamma_mom_neutral, r.Gamma_molecule, ...12
 r.core2sol_flux, r.core2sol_mol, r.sol2extern_ion_flux,  ...
 r.Source_n, r.Source_v, r.Source_Q, r.Source_neutral, r.Source_vn, r.Source_molecule,  ...21
 r.extern_neutral_density,  ... 22
 r.extern_neutral_flux, r.sol2extern_flux, r.extern2sol_flux, r.tar2extern_flux,  ... 26
 r.extern2core_flux, r.sum_sol2extern_ion_flux, ... 28
 r.Source_extern, r.neutral_pump, ...30
 r.extern_molecule_density,  ...31
 r.extern_molecule_flux, r.sol2extern_mol, r.extern2sol_mol, r.tar2extern_mol, ...35
 r.extern2core_mol, r.sum_sol2extern_ion_mol,  ...37
 r.Source_extern_mol, r.molecule_pump, ...39
 r.core_density, r.Gamma_core2sol, r.sol2core_flux, r.sol2core_mol, r.Source_core, ...
 r.end_time, r.core_neutral_density]...
 = clib.libdiv1d.run_div1d_(	...
 a.density, a.velocity, a.temperature, a.neutral_density, a.neutral_velocity, a.molecule, ... 6
 r.Gamma_n, r.Gamma_mom, r.q_parallel, r.Gamma_neutral, r.Gamma_mom_neutral, r.Gamma_molecule, ...12
 r.core2sol_flux, r.core2sol_mol, r.sol2extern_ion_flux,  ...
 r.Source_n, r.Source_v, r.Source_Q, r.Source_neutral, r.Source_vn, r.Source_molecule,  ...21
 r.extern_neutral_density,  ... 22
 r.extern_neutral_flux, r.sol2extern_flux, r.extern2sol_flux, r.tar2extern_flux,  ... 26
 r.extern2core_flux, r.sum_sol2extern_ion_flux, ... 28
 r.Source_extern, r.neutral_pump, ...30
 r.extern_molecule_density,  ...31
 r.extern_molecule_flux, r.sol2extern_mol, r.extern2sol_mol, r.tar2extern_mol, ...35
 r.extern2core_mol, r.sum_sol2extern_ion_mol,  ...37
 r.Source_extern_mol, r.molecule_pump, ...39
 r.core_density, r.Gamma_core2sol, r.sol2core_flux, r.sol2core_mol, r.Source_core, ...44
 r.start_time, r.end_time, r.nout, r.delta_t,  ...48
 r.imp_con, r.neu, r.dneu, r.nb, r.mb, ... 53 
 r.gas, r.rec, r.qpar_x, r.red_frc, ... 57
 r.Q_core, r.Gamma_core, r.core_neutral_density , ...60
 r.neutral_puff, r.molecule_puff, r.core_fuelling); ...63

r.time = r.start_time:r.delta_t*r.nout:r.end_time-r.delta_t*r.nout;
r.X = a.x;
r.Xcb = a.xcb;

%% r = div1d_lib_runvars(i.intinnum(1));
addpath(genpath('/home/unix/derks/Desktop/projects/dynamics/toolbox'))
tmp{1}.divout = r;
plotdiv1d_v600(tmp{1}.divout,input_struct)
for i_step = 1:10
    tmp{i_step+1}.divout = div1d_runner(tmp{i_step}.divout); % step 
    plotdiv1d_v600(tmp{i_step+1}.divout,input_struct)
end
plotdiv1d_v600(tmp{1}.divout,input_struct)
%
disp('M finished calling div1d library, got return values')

%% change grid and profiles on the go:
[a.density, a.velocity, a.temperature, a.neutral_density, a.neutral_velocity, a.molecule, ...
		           a.x, a.xcb, a.delta_x, a.delta_xcb, a.b_fieldd, a.b_field_cb, a.b_trans, a.b_trans_cb, ...
			   a.r_cc, a.r_cb, a.area_extern, a.sintheta_cc, a.sintheta_cb, a.sol_width_pol, a.sol_width_pol_cb, a.volumes,...
			   a.gas_puff_profile, a.core_source_profile_q, a.core_source_profile_n,...
		       a.i_omp, a.i_xpoint, a.i_baffle, a.mid_point, a.x_omp] ...
    =clib.libdiv1d.initialize_div1d_arrays_(a.density, a.velocity, a.temperature, a.neutral_density, a.neutral_velocity, a.molecule, ...
		           a.x, a.xcb, a.delta_x, a.delta_xcb, a.b_field, a.b_field_cb, a.b_trans, a.b_trans_cb, ...
			   a.r_cc, a.r_cb, a.area_extern, a.sintheta_cc, a.sintheta_cb, a.sol_width_pol, a.sol_width_pol_cb, a.volumes,...
			   a.gas_puff_profile, a.core_source_profile_q, a.core_source_profile_n,...
		           a.init_grid_fortran, a.init_prof_fortran, i.intinnum(1), a.i_omp, a.i_xpoint, a.i_baffle, a.mid_point, a.x_omp);
disp('M finished re-setting the grid (this makes the grid variable! ')