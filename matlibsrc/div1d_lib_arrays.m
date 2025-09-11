function [o] = div1d_lib_arrays(nx)
o = struct;
o.density( nx) = 1.0;
o.velocity( nx) = 1.0;
o.temperature( nx) = 1.0;
o.neutral_density( nx) = 1.0; 
o.neutral_velocity( nx) = 1.0;
o.molecule( nx) = 1.0;
o.x( nx) = 1.0;
o.xcb( nx+1) = 1.0;
o.delta_x( nx-1) = 1.0;
o.delta_xcb( nx) = 1.0;
o.b_field( nx) = 1.0;
o.b_field_cb( nx+1) = 1.0;
o.b_trans( nx) = 1.0;
o.b_trans_cb( nx+1) = 1.0;
o.r_cc( nx) = 1.0;
o.r_cb( nx+1) =1.0;
o.area_extern( nx) = 1.0;
o.sintheta_cc( nx) = 1.0;
o.sintheta_cb( nx+1) = 1.0;
o.sol_width_pol( nx) = 1.0;
o.sol_width_pol_cb( nx+1) =1.0;
o.volumes(nx) = 1.0;
o.gas_puff_profile(nx) = 1.0;
o.core_source_profile_q(nx) = 1.0;
o.core_source_profile_n(nx) = 1.0;
o.init_grid_fortran = 1;
o.init_prof_fortran = 1;
o.nx = 200;
o.i_omp = 1;
o.i_xpoint(2) = 1;
o.i_baffle(2) = 1;
o.mid_point = 1.0;
o.x_omp = 1.0;


end