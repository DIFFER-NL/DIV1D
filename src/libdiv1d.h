#ifndef libdiv1d_h__
#define libdiv1d_h__

#include <stdio.h> // namespace with i/o stream
#include <math.h>
#include <stdlib.h>
// * requests adress of variable to pass by reference. 
// https://stackoverflow.com/questions/11428526/passing-a-matrix-in-a-function-c 



// #ifdef INITIALIZE_SETTINGS
extern void initialize_div1d_settings_(double *floatinnum, int *intinnum, signed char *loginnum, double *floatinphys, int *intinphys, int *call_from_extern);
// #ifdef INITIALIZE_DIV1D 
extern void initialize_div1d_arrays_(double *E_density, double *E_velocity, double *E_temperature, double *E_neutral, double *E_neutral_velocity, double *E_molecule,
		           double *E_x, double *E_xcb, double *E_delta_x, double *E_delta_xcb, double *E_B_field, double *E_B_field_cb, double *E_B_trans, double *E_B_trans_cb,
			   double  *E_R_cc, double *E_R_cb, double *E_Area_extern, double *E_sintheta_cc, double *E_sintheta_cb, double *E_sol_width_pol, double *E_sol_width_pol_cb, double *E_volumes, 
			   double *E_gas_puff_profile, double *E_core_source_profile_q, double *E_core_source_profile_n,
		           int *init_grid_fortran, int *init_prof_fortran, int *Nx, int *E_i_omp, int *E_i_Xpoint, int *E_i_baffle, double *E_mid_point, double *E_X_omp);	 
// #ifdef RUN_DIV1D
// init_div1d_matlab ( var sizes, Nx )
//	int nx
//	malloc denisty (NX) 

extern void run_div1d_(	double *E_density, double *E_velocity, double *E_temperature, double *E_neutral, double *E_neutral_velocity, double *E_molecule, // 6
		       	double *E_Gamma_n, double *E_Gamma_mom, double *E_q_parallel, double *E_Gamma_neutral, double *E_Gamma_mom_neutral, double *E_Gamma_molecule, //12
			double *E_core2sol_flux, double *E_core2sol_mol, double *E_sol2extern_ion_flux,  //15
			double *E_Source_n, double *E_Source_v, double *E_Source_Q, double *E_Source_neutral, double *E_Source_vn, double *E_Source_molecule,  //21
			double *E_extern_neutral_density,  // 22
			double *E_extern_neutral_flux, double *E_sol2extern_flux, double *E_extern2sol_flux, double *E_tar2extern_flux,  // 26
			double *E_extern2core_flux, double *E_sum_sol2extern_ion_flux, // 28
			double *E_Source_extern, double *E_neutral_pump, //30
			double *E_extern_molecule_density,  //31
			double *E_extern_molecule_flux, double *E_sol2extern_mol, double *E_extern2sol_mol, double *E_tar2extern_mol, //35
			double *E_extern2core_mol, double *E_sum_sol2extern_ion_mol,  //37
			double *E_Source_extern_mol, double *E_molecule_pump, //39
			double *E_core_density, double *E_Gamma_core2sol, double *E_sol2core_flux, double *E_sol2core_mol, double *E_Source_core, //44
		        double *start_time, double *end_time, int *nout, double *delta_t,  //48
		        double *E_imp_con, double *E_neu, double *E_dneu, double *E_nb, double *E_mb, // 53
		        double *E_gas, double *E_rec, double *E_qpar_x, double *E_red_frc, // 57
		        double *E_Q_core, double *E_Gamma_core, double *E_core_neutral_density , //60
		        double *E_neutral_puff, double *E_molecule_puff, double *E_core_fuelling); //63
//#endif // RUN_DIV1D
//#endif // INITIALIZE_DIV1D
//#endif // INITIALIZE_SETTINGS

#endif // libdiv1d_h__
