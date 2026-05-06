% build the C interface to the DIV1D library

% run this script from the div1d folder
clear vars;
addpath(genpath('/lustre/isaac24/proj/UTK0325/leekuanwei/div1d_astra'))
addpath(genpath('/lustre/isaac24/proj/UTK0325/leekuanwei/div1d_astra/div1d-matlab-wrapper'))
addpath(genpath('/lustre/isaac24/proj/UTK0325/leekuanwei/div1d_astra/div1d_runs/solps-1d-toolbox'))
addpath(genpath('/lustre/isaac24/proj/UTK0325/leekuanwei/div1d_astra/div1d_runs/solps-iter_processing'))
system('./build_libdiv1d.sh');
filePath = pwd;
addpath('obj');
addpath('src');
addpath('matlibsrc')
addpath('namelist')

%run C11_get_MAST_solps_information.m
%close all

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
libpath = '/lustre/isaac24/proj/UTK0325/leekuanwei/div1d_astra';
libDef = definelibdiv1d();
% Rsummary(definelibdiv1d(libpath))
summary(libDef)
libDef = definelibdiv1d_extend_test(libDef);
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

%% Core Particle balance
simpath = '/lustre/isaac24/proj/UTK0325/leekuanwei/div1d_astra/div1d_runs/cntrl_testcase/div1d_output.txt';
[output,input] = div1dread_v600(simpath);

%plotdiv1d_v600(output,input)
o = get_last_output(output);
g = input.grid;
g.delta_x = diff(g.x);
g.delta_xcb = diff(g.xcb);

% ato_pump = [0 0 5e2 5e2 5e2];
% mol_pump = [0 0 1e2 1e2 1e2];
% %ato_pump = [0 0 0 0 0];
% %mol_pump = [0 0 0 0 0];
% core_fueling = 1000;
% wall_ass_prob = 0.01;
% [sim, set] = get_div1d_chamber_params(input,output,'ato_pump',ato_pump,'mol_pump',mol_pump,...
%     'tau_core',7e-2,'core_fuelling',core_fueling,'wall_ass_prob',wall_ass_prob);
%% define DIV1D object 
% settings
%addpath('/home/unix/derks/Desktop/codes/div1d_library/libdiv1d')

%i = div1d_lib_inputs;
i = div1d_lib_inputs_test(input,o);
i.call_from_extern = 1;
%
if init_div1d_settings ==0
clib.libdiv1d.initialize_div1d_settings_( ...
            i.floatinnum, i.intinnum, i.loginnum, ...
            i.floatinphys,i.intinphys, i.call_from_extern  );
init_div1d_settings = 0
 else
 disp('M cannot initialize div1d twice, segmentation fault for second allocations')
 end

% works! 
%%
if init_div1d_arrays == 0
a = div1d_lib_arrays(i.intinnum(1)); 
a.init_grid_fortran = 0;
a.init_prof_fortran = 0;
% [a.density, a.velocity, a.temperature, a.neutral_density, a.neutral_velocity, a.molecule, ...
% 		           a.x, a.xcb, a.delta_x, a.delta_xcb, a.b_fieldd, a.b_field_cb, a.b_trans, a.b_trans_cb, ...
% 			   a.r_cc, a.r_cb, a.area_extern, a.sintheta_cc, a.sintheta_cb, a.sol_width_pol, a.sol_width_pol_cb, a.volumes,...
% 			   a.gas_puff_profile, a.core_source_profile_q, a.core_source_profile_n,...
% 		       a.i_omp, a.i_xpoint, a.i_baffle, a.mid_point, a.x_omp, a.A_int, a.A_wet] ...
%     =clib.libdiv1d.initialize_div1d_arrays_(a.density, a.velocity, a.temperature, a.neutral_density, a.neutral_velocity, a.molecule, ...
% 		           a.x, a.xcb, a.delta_x, a.delta_xcb, a.b_field, a.b_field_cb, a.b_trans, a.b_trans_cb, ...
% 			   a.r_cc, a.r_cb, a.area_extern, a.sintheta_cc, a.sintheta_cb, a.sol_width_pol, a.sol_width_pol_cb, a.volumes,...
% 			   a.gas_puff_profile, a.core_source_profile_q, a.core_source_profile_n,...
% 		           a.init_grid_fortran, a.init_prof_fortran, i.intinnum(1), a.i_omp, a.i_xpoint, a.i_baffle, a.mid_point, a.x_omp, g.a_int, g.a_wet);

[a.density, a.velocity, a.temperature, a.neutral_density, a.neutral_velocity, a.molecule, ...
		           a.x, a.xcb, a.delta_x, a.delta_xcb, a.b_field, a.b_field_cb, a.b_trans, a.b_trans_cb, ...
			   a.r_cc, a.r_cb, a.area_extern, a.sintheta_cc, a.sintheta_cb, a.sol_width_pol, a.sol_width_pol_cb, a.volumes,...
			   a.gas_puff_profile, a.core_source_profile_q, a.core_source_profile_n,...
		       a.i_omp, a.i_xpoint, a.i_baffle, a.mid_point, a.x_omp, a.a_int, a.a_wet] ...
    =clib.libdiv1d.initialize_div1d_arrays_(o.density, o.velocity, o.temperature, o.neutral_density, o.neutral_velocity, o.molecule, ...
		           g.x, g.xcb, g.delta_x, g.delta_xcb, g.b_field, g.b_field_cb, g.b_trans, g.b_trans_cb, ...
			   g.r_cc, g.r_cb, g.area_extern, g.sintheta_cc, g.sintheta_cb, g.sol_width_pol, g.sol_width_pol_cb, g.volumes,...
			   g.gas_puff_profile, g.e_core_source_profile_q, g.e_core_source_profile_n,...
		           a.init_grid_fortran, a.init_prof_fortran, i.intinnum(1), g.i_omp, g.i_xpoint, g.i_baffle, g.mid_point, g.x_omp, g.a_int, g.a_wet);
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
%a = div1d_lib_arrays(i.intinnum(1)); 
%% now all has been initialized we can get the div1d input struct from the
% echo file.
% echo_path = '/lustre/isaac24/proj/UTK0325/leekuanwei/div1d_astra/echo_div1d_inputs.txt';
% input_struct = read_echo_div1d_inputs(echo_path);
% input_struct.grid = a;
%% Core Particle balance
% ato_pump = [0 0 5e2 5e2 5e2];
% mol_pump = [0 0 1e2 1e2 1e2];
% %ato_pump = [0 0 0 0 0];
% %mol_pump = [0 0 0 0 0];
% core_fueling = 1000;
% wall_ass_prob = 0.01;
% [sim, set] = get_div1d_chamber_params(input_struct,r,'ato_pump',ato_pump,'mol_pump',mol_pump,...
%     'tau_core',7e-2,'core_fuelling',core_fueling,'wall_ass_prob',wall_ass_prob);
%%

r = div1d_runstruct(a.nx,i.intinnum(2), i);
r.density = o.density;
r.velocity = o.velocity;
r.neutral_density = o.neutral_density;
r.neutral_velocity = o.neutral_velocity;
r.molecule = o.molecule;
r.temperature = o.temperature;
r.core_density = o.core_density;

r.Q_core = 1 * input.physics.q_core * ones(1,r.nout);
r.Gamma_core = input.physics.gamma_core * ones(1,r.nout);
r.core_fuelling = input.physics.core_fuelling * ones(1,r.nout);
r.neutral_pump = input.physics.pump_rate_n;
r.molecule_pump = input.physics.pump_rate_m;
%r.neutral_puff = zeros(5,10);
%r.molecule_puff = zeros(5,10);
%r.neutral_pump = [0 0 0 0 0];
%r.molecule_pump = [0 0 0 0 0];
r.delta_t = 1e-5;
%r.core_fuelling = 1000 * ones(1,r.nout);
r.rec = input.physics.recycling * ones(1,r.nout);
r.qpar_x = input.physics.q_parx * ones(1,r.nout);

%r.imp_con = repmat([0.002 ,0,0, 0, 0 ],5,2);

r.molecule_puff = input.physics.puff_rate_molecule' * ones(1,r.nout);
r.neutral_puff = input.physics.puff_rate_neutral' * ones(1,r.nout);

% r.nb(1,:) = o.extern_neutral_density(1);
% r.nb(2,:) = o.extern_neutral_density(2);
% r.nb(3,:) = o.extern_neutral_density(3);
% r.nb(4,:) = o.extern_neutral_density(4);
% r.nb(5,:) = o.extern_neutral_density(5);
% 
% r.mb(1,:) = o.extern_molecule_density(1);
% r.mb(2,:) = o.extern_molecule_density(2);
% r.mb(3,:) = o.extern_molecule_density(3);
% r.mb(4,:) = o.extern_molecule_density(4);
% r.mb(5,:) = o.extern_molecule_density(5);
%r.molecule_puff = input.physics.puff_rate_molecule' * ones(1,r.nout);
%r.neutral_puff = input.physics.puff_rate_neutral' * ones(1,r.nout);

% [r.density, r.velocity, r.temperature, r.neutral_density, r.neutral_velocity, r.molecule, ... 6
%  r.Gamma_n, r.Gamma_mom, r.q_parallel, r.Gamma_neutral, r.Gamma_mom_neutral, r.Gamma_molecule, ...12
%  r.core2sol_flux, r.core2sol_mol, r.sol2extern_ion_flux,  ...
%  r.Source_n, r.Source_v, r.Source_Q, r.Source_neutral, r.Source_vn, r.Source_molecule,  ...21
%  r.extern_neutral_density,  ... 22
%  r.extern_neutral_flux, r.sol2extern_flux, r.extern2sol_flux, r.tar2extern_flux,  ... 26
%  r.extern2core_flux, r.sum_sol2extern_ion_flux, ... 28
%  r.Source_extern, r.neutral_pump, ...30
%  r.extern_molecule_density,  ...31
%  r.extern_molecule_flux, r.sol2extern_mol, r.extern2sol_mol, r.tar2extern_mol, ...35
%  r.extern2core_mol, r.sum_sol2extern_ion_mol,  ...37
%  r.Source_extern_mol, r.molecule_pump, ...39
%  r.core_density, r.Gamma_core2sol, r.sol2core_flux, r.sol2core_mol, r.Source_core, ...
%  r.end_time, r.core_neutral_density]...
%  = clib.libdiv1d.run_div1d_(	...
%  o.density, o.velocity, o.temperature, o.neutral_density, o.neutral_velocity, o.molecule, ... 6
%  o.Gamma_n, o.Gamma_mom, o.q_parallel, o.Gamma_neutral, o.Gamma_mom_neutral, o.Gamma_molecule, ...12
%  o.core2sol_flux, o.core2sol_mol, r.sol2extern_ion_flux,  ...
%  o.Source_n, o.Source_v, o.Source_Q, o.Source_neutral, o.Source_vn, o.Source_molecule,  ...21
%  o.extern_neutral_density,  ... 22
%  o.extern_neutral_flux, o.sol2extern_flux, o.extern2sol_flux, o.tar2extern_flux,  ... 26
%  o.extern2core_flux, o.sum_sol2extern_ion_flux, ... 28
%  r.Source_extern, r.neutral_pump, ...30
%  o.extern_molecule_density,  ...31
%  o.extern_molecule_flux, o.sol2extern_mol, o.extern2sol_mol, o.tar2extern_mol, ...35
%  o.extern2core_mol, o.sum_sol2extern_ion_mol,  ...37
%  r.Source_extern_mol, r.molecule_pump, ...39
%  o.core_density, o.Gamma_core2sol, o.sol2core_flux, o.sol2core_mol, o.Source_core, ...44
%  r.start_time, r.end_time, r.nout, r.delta_t,  ...48
%  r.imp_con, r.neu, r.dneu, r.nb, r.mb, ... 53 
%  r.gas, r.rec, r.qpar_x, r.red_frc, ... 57
%  r.Q_core, r.Gamma_core, r.core_neutral_density , ...60
%  r.neutral_puff, r.molecule_puff, r.core_fuelling); ...63
r.mb = r.mb';
r.nb = r.nb';
r.imp_con = r.imp_con';
r.neutral_puff = r.neutral_puff';
r.molecule_puff = r.molecule_puff';

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
 r.density, r.velocity, r.temperature, r.neutral_density, r.neutral_velocity, r.molecule, ... 6
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


%% Core Particle balance
% ato_pump = [0 0 5e2 5e2 5e2];
% mol_pump = [0 0 1e2 1e2 1e2];
% %ato_pump = [0 0 0 0 0];
% %mol_pump = [0 0 0 0 0];
% core_fueling = 1000;
% wall_ass_prob = 0.01;
% [sim, set] = get_div1d_chamber_params(input_struct,r,'ato_pump',ato_pump,'mol_pump',mol_pump,...
%     'tau_core',7e-2,'core_fuelling',core_fueling,'wall_ass_prob',wall_ass_prob);
%% r = div1d_lib_runvars(i.intinnum(1));
% rundir = '/fusion/projects/solps-results/leekuanwei/runs/MASTU_cd_Hmode_DC/mastu_Hmode_1e22_mid/'; 
% load([rundir, '/solpsdata_gijs.mat'])
% 
% for kk = 1:40
%     for kj = 1:148
% % volumes1 = 
%         volumes(kj,kk) = areasquare(shot.geom.r(kj,kk,:),shot.geom.z(kj,kk,:))*shot.geom.cr(kj,kk)*2*pi;
%     end
% end
% 
% shot.geom.volumes = volumes;
% shot.geom.sep = 22; 
% iXp = [[shot.geom.leftcut'+[0 2]] [flip(shot.geom.rightcut)'+[1 2]]];
% shot.geom.omp = 63 + iXp(2); 
% dev = 1;
% tmp2 = get_solps_profiles(shot,'dev',dev,'cases',{'maksi'},'wait',0, 'iXp', [12,44,87,118]); 
% solps = tmp2.maksi;
%addpath(genpath('/home/unix/derks/Desktop/projects/dynamics/toolbox'))
%r = tmp{end}.divout;
%tmp = {};
%%
tmp = {};
tmp{1}.divout = r;
for i_step = 1:1
    tmp{i_step+1}.divout = div1d_runner(tmp{i_step}.divout); % step 
    %te_hist(i) = tmp{i_step+1}.divout.temperature(end);
    %plotdiv1d_v600(tmp{i_step+1}.divout,input_struct)
end
r = tmp{end}.divout;

%%
% for i = 1:11
%     te_hist(i) = tmp{i}.divout.temperature(end);
% 
% end   
% 
% plot(te_hist)
%%
figure(1)
plot(r.X,r.temperature,Color='blue')
hold on
plot(o.X,o.temperature,Color='red')
% figure(2)
% plot(r.Xcb,r.q_parallel,Color='blue')
% hold on
% plot(o.Xcb,o.q_parallel,Color='red')
% plot(r.Xcb,r.q_parallel .* g.b_field_cb,Color= 'green')
%%

params.Kp_Te   = 5e20;
params.Ki_Te   = 1.0;

params.Kp_qpar = 0.5;
params.Ki_qpar = 0.1;

params.Kp_jsat = 0;

params.Kp_Te_det   = 10.0;
params.Ki_Te_det   = 2.0;
params.Kp_qpar_det = 0.2;

params.jsat_rollover_thresh = 1e4;   % A/m^2/s threshold
params.Imax = 10.0;

params.gas_puff_min = 0.0;
params.gas_puff_max = 4e22;
params.dgas_max     = 5.0e20;

params.gas_puff_init = 0;
tmp = {};
tmp{1}.divout = r;

%plotdiv1d_v600(tmp{1}.divout,input_struct)
Te_ref = 5;
qpar_ref = 2e6;
jsat_ref = 1e4;
ctrl = {};
dt = r.delta_t;
gas_background = r.molecule_puff(1,4);
Te = r.temperature(end);
qpar = r.q_parallel(end);
gas = r.molecule_puff(1,4)-gas_background;
jsat = r.Gamma_n(end) * 1.6-19;
Te_hist(1) = Te;
qpar_hist(1) = qpar;
gas_hist(1) = gas;
jsat_hist(1) = jsat;
u_hist(1) = 0;
Time(1) = 0;

for k = 1:45

    [gas_puff, ctrl,u] = detachment_pid_controller( ...
        Te_ref, qpar_ref, jsat_ref, ...
        Te, qpar, jsat, ...
        ctrl, params, dt);


    r.molecule_puff(:,4) = gas_background * ones(r.nout,1) + gas_puff*ones(r.nout,1);
    tmp = {};
    tmp{1}.divout = r;
    for i_step = 1:1
        tmp{i_step+1}.divout = div1d_runner(tmp{i_step}.divout); % step 
        %plotdiv1d_v600(tmp{i_step+1}.divout,input_struct)
    end
    r = tmp{end}.divout;
    Te = r.temperature(end);
    qpar = r.q_parallel(end);
    gas = r.molecule_puff(1,4)-gas_background;
    jsat = r.Gamma_n(end) * 1.6-19;
    Te_hist(k+1) = Te;
    qpar_hist(k+1) = qpar;
    gas_hist(k+1) = gas;
    jsat_hist(k+1) = jsat;
    u_hist(k+1) = u;

    Time(k+1) = r.delta_t * r.nout * k;
end
ctrlvars = {};
ctrlvars.Time = Time;
ctrlvars.Te = Te_hist;
ctrlvars.qpar = qpar_hist;
ctrlvars.gas_puff = gas_hist;
ctrlvars.jsat = jsat_hist;
ctrlvars.u = u_hist;

save('ctrlvars.mat',"ctrlvars")
%r = tmp{end}.divout;
% %plotdiv1d_v600(tmp{1}.divout,input_struct)
%
disp('M finished calling div1d library, got return values')

%%
% for i = 1:5
%     core_density(i) = tmp{i}.divout.core_density;
%     background_density(1,i) = tmp{i}.divout.extern_neutral_density(3);
%     background_density(2,i) = tmp{i}.divout.extern_neutral_density(4);
%     background_density(3,i) = tmp{i}.divout.extern_neutral_density(5);
% 
%     background_density_mol(1,i) = tmp{i}.divout.extern_molecule_density(3);
%     background_density_mol(2,i) = tmp{i}.divout.extern_molecule_density(4);
%     background_density_mol(3,i) = tmp{i}.divout.extern_molecule_density(5);
%     temperature(i) = tmp{i}.divout.temperature(end); 
% end
% figure(1)
% plot(core_density)
% 
% figure(2)
% plot(background_density(1,:))
% 
% figure(3)
% plot(background_density_mol(1,:))
% 
% figure(4)
% plot(temperature)
% for i = 1:10
%     plot(i,tmp{i}.divout.)
% end    

%% Post process
% r.sum_extern2core_flux = sum(r.extern2core_flux);
% r.sum_extern2core_mol = sum(r.extern2core_mol);
% [out] = process_div1d_output_v600(r,input,'plot',1,'hold',1);
% 
% plotdiv1d_v600(output,input)
% hold on
% plotdiv1d_v600(r,input)
% plotdiv1d_v600(tmp{11}.divout,input)

%% Read TXT File
% importeddata = readmatrix('density_last.txt');
% importeddata = readmatrix('Source_n_first.txt');
% density = [];
% for i = 1:length(importeddata(:,1))
%     density = [density importeddata(i,:)];
% end
% density(1) = [];
% figure(1)
% plot(o.Source_n)
% hold on
% plot(density)
% legend('standalone','matlab')
% 
% figure(2)
% plot(output.density(9,:))
% hold on
% plot(o.density)
% legend('Second last','last')
%% Test plots
% figure(1)
% plot(o.Gamma_neutral,LineWidth=3, Color='b')
% hold on
% plot(r.Gamma_neutral,LineWidth=3,Color='r')
% %plot(r.q_parallel,LineWidth=3,Color='green')
% %hold on
% % figure(2)
% % plot(r.Gamma_mom .* g.b_field_cb,LineWidth=3,Color='green')
% % 
%legend('standalone','matlab')
%%
%Read Source_n.txt and separate data blocks by the string "Source_n:"
% filename = 'fort.210';
% filename_standalone = '/home/leekuanwei/div1d_astra/div1d_runs/mid_2e22/fort.210';
% % 
%  variable = 'Source_Q';
% % 
% data_block = readtxt_file(filename,variable);
% data_block_standalone = readtxt_file(filename_standalone,variable);
% % 
% figure(10)
% 
% plot(data_block{1},LineWidth=3,Color='b')
% hold on
% plot(data_block_standalone{end},LineWidth=3,Color='r')
% %plot(o.Source_neutral(end,:))
% legend('Matlab','Standalone')

% figure(10)
% 
% plot(data_block{end},LineWidth=3,Color='b')
% hold on
% plot(o.Source_n,LineWidth=3,Color='r')
% legend('Matlab','Standalone')

%%
%run C11_get_MAST_solps_information.m
%close all
% titlestr = 'MAST-U';
% L = input.physics.l;
% r.q_parallel = r.q_parallel .* g.b_field_cb;
% xticks = round(linspace(0,L,10),0);
% plot_div1d_profiles(r,input,'solps',solps,'figtight',1,'title',titlestr,'xlimits',[0 L],...
%                 'xp',L*1.16,'generalposition',[100 100 370 740],'tf_form',[0.2 0.04 0.06 0.04],'pbline',1,...
%                 'xticks',xticks,'Xpoint',input.grid.i_xpoint(2),'fignum',1,'hold',1,'plotimpurities',0,...
%                 'ytic_a',[10^15 10^16 10^17 10^18],'ylim_a', [0.6*10^16 9*10^18],'ylim_v',[0 25],'ytic_v',[0 10 20 30],...
%                 'ytic_t',[0 20 40 60],'ylim_t',[0 50],'ylim_n',[3*10^18 6e19],...
%                 'ytic_m',[10^15 10^16 10^17 10^18],'ylim_m',[5*10^14 5*10^19],'ylim_vn', [-1,15], ...
%                 'plotbackground',1,'pxline',1, 'FontSize', 13, 'plotmolecules', 1, 'plotvn', 1);


%% change grid and profiles on the go:
% [a.density, a.velocity, a.temperature, a.neutral_density, a.neutral_velocity, a.molecule, ...
% 		           a.x, a.xcb, a.delta_x, a.delta_xcb, a.b_fieldd, a.b_field_cb, a.b_trans, a.b_trans_cb, ...
% 			   a.r_cc, a.r_cb, a.area_extern, a.sintheta_cc, a.sintheta_cb, a.sol_width_pol, a.sol_width_pol_cb, a.volumes,...
% 			   a.gas_puff_profile, a.core_source_profile_q, a.core_source_profile_n,...
% 		       a.i_omp, a.i_xpoint, a.i_baffle, a.mid_point, a.x_omp] ...
%     =clib.libdiv1d.initialize_div1d_arrays_(a.density, a.velocity, a.temperature, a.neutral_density, a.neutral_velocity, a.molecule, ...
% 		           a.x, a.xcb, a.delta_x, a.delta_xcb, a.b_field, a.b_field_cb, a.b_trans, a.b_trans_cb, ...
% 			   a.r_cc, a.r_cb, a.area_extern, a.sintheta_cc, a.sintheta_cb, a.sol_width_pol, a.sol_width_pol_cb, a.volumes,...
% 			   a.gas_puff_profile, a.core_source_profile_q, a.core_source_profile_n,...
% 		           a.init_grid_fortran, a.init_prof_fortran, i.intinnum(1), a.i_omp, a.i_xpoint, a.i_baffle, a.mid_point, a.x_omp);
% disp('M finished re-setting the grid (this makes the grid variable! ')

% simpath = '/fusion/projects/solps-results/leekuanwei/runs/div1d_astra/div1d_runs/mid_1e22_6e22_mol_div/div1d_output.txt';
% [output,input] = div1dread_v600(simpath);
% plotdiv1d_v600(output,input)