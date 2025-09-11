// Example of a C script calling DIV1D
//
// everything is passed by reference and should be a pointer
// single values are made pointer by initializing int *Nx = 300;
// arrays are always initialized as lists of pointers double array[30];




#include <stdio.h> // namespace with i/o stream
#include <math.h>
#include <stdlib.h>

#define INITIALIZE_SETTINGS
#define INITIALIZE_DIV1D
#define RUN_DIV1D

#ifdef INITIALIZE_SETTINGS
#include <libdiv1d.h>
#endif

int main(void)
{
 int i;
 // settings
 #ifdef INITIALIZE_SETTINGS
 
 // numerics settings
 double floatinnum[90];
for(i=0;i<90;++i) 
{
floatinnum[i] = 0.0; 
}
 int intinnum[30];
for(i=0;i<30;++i) 
{
intinnum[i] = 0; 
}

signed char loginnum[10] = {0,0,1,0,0,0,0,0,0,0};

double floatinphys[141];
for(int i=0;i<141; ++i) 
{
floatinphys[i] = 0.0; 
}
int intinphys[31];
for(int i=0;i<31; ++i) 
{
intinphys[i] = 0; 
}

// set values
intinnum[0] = 100; // Nx
intinnum[1] = 10; // nout
intinnum[2] = 227; // method
intinnum[3] = 10000; // itate_mod
intinnum[4] = 10000; // max_step
intinnum[5] = 100; // max attempts
intinnum[6] = -100; // nzswag
intinnum[7] = 1; // evolve density
intinnum[8] = 1; // evolve momentum
intinnum[9] = 1; // evolve energy
intinnum[10] = 1; // evolve neutral
intinnum[11] = 1; // evolve neutral momentum
intinnum[12] = 1; // evolve molecule
intinnum[13] = 1; // evolve background
intinnum[14] = 1; // evolve core
intinnum[15] = 0; // evolve core neutral
intinnum[16] = 1; // mol dens model
intinnum[17] = 1; // D_new harmonic average of temperature for diffusion
// norms
floatinnum[0] = 1.0*pow(10,19); // density norm
floatinnum[1] = 1.0; // temperature norm
floatinnum[2] = 1.0*pow(10,4); // velocity norm
floatinnum[3] = 1.0; // neutral norm
// equation switch
floatinnum[9] = 1.0; //  density source
floatinnum[10] = 1.0; // momentum source
floatinnum[11] = 1.0; // energy source
floatinnum[12] = 1.0; // neutral source
floatinnum[13] = 1.0; // neutral momentum source
floatinnum[14] = 1.0; // molecule source
// switches energy terms
floatinnum[19] = 1.0; // conv heat
floatinnum[20] = 1.0; // energy flux
floatinnum[21] = 1.0; // energy compr
floatinnum[22] = 0.0; // only core source Q
// switches atomic reaction rates + impurity
floatinnum[29] = 1.0; // cx
floatinnum[30] = 1.0; // rec
floatinnum[31] = 1.0; // rec ene
floatinnum[32] = 1.0; // ion 
floatinnum[33] = 1.0; // exc 
floatinnum[34] = 1.0; // impurity  
floatinnum[35] = 1.0; // mom atoms
floatinnum[36] = 1.0; // mom molecules
// switches molecular rates
floatinnum[39] = 1.0; //  diss mol
floatinnum[40] = 1.0; //  ion mol
floatinnum[41] = 1.0; //  cx mol
floatinnum[42] = 1.0; //  diss h2pl 
floatinnum[43] = 1.0; //  dis rec h2pl
floatinnum[44] = 1.0; //  dis ion h2pl
floatinnum[45] = 1.0; //  energy molecules
floatinnum[46] = 1.0; //  diss att
floatinnum[47] = 1.0; //  cx hmin
floatinnum[48] = 1.0; // ion hmin
// numerical settings
floatinnum[59] = 0.05; // dxmin
floatinnum[60] =pow(10,-5); //  delta t
floatinnum[61] = pow(10,-6); //  abstol
floatinnum[62] = pow(10,-6); //  reltol
floatinnum[63] = 5; // vis
floatinnum[64] = 0.5; // cent diff
floatinnum[65] = 1.0; // lax switch


 // physics settings
// phys integers
intinphys[0] = 5  ;// num imp
intinphys[1] = 6  ;// imp1
intinphys[2] = 7 ;// imp2
intinphys[3] = 0 ;// imp3
intinphys[4] = 0 ;// imp4
intinphys[5] = 0 ;// imp5
intinphys[6] = 0  ;// switch imp dist
intinphys[7] = 0 ;// widt pfr

 // phys floats
// geometry, grid distributions
floatinphys[0] = 20.0 ;// L
floatinphys[1] = 10.0 ;// L cor sol
floatinphys[2] = 0.0 ;// x cor sol
floatinphys[3] = 0.1 ;// sintheta
floatinphys[4] = 90.0 ;// pol tar ang 1
floatinphys[5] = 90.0 ;// pol tar ang 2
floatinphys[6] = 1.1*pow(10,-2) ;// sol wid omp 
floatinphys[7] = 0.0 ;// loc omp 
floatinphys[8] = 1.4 ;// maj rad

floatinphys[9] = 1.0 ;// alpha core prof q
floatinphys[10] = 1.0*pow(10,-3) ;// alpha core prof n
floatinphys[11] = 2.0 ;// flux expansion
floatinphys[12] = 3.0 ;// trans expansion
floatinphys[13] = 0.0 ;// gas puf loc
floatinphys[14] = pow(10,20) ;// gas puf wid
floatinphys[15] = 0.4500 ;// sig nb backgrounds
floatinphys[16] = 0.0; // L_baffle

// sol influx settings
floatinphys[19] = 0.0 ;// gamma_x
floatinphys[20] = 0.0 ;// qparx
floatinphys[21] = pow(10,20) ;// gamma core
floatinphys[22] = 2.2*pow(10,5) ;// q core
floatinphys[23] = 0.0 ;// density ramp rate
floatinphys[24] = 0.0 ;// gas puf src

// Q_CORE =  2.21735701291665e+05
//  GAMMA_CORE = 1.0000e+20

// initial values
floatinphys[29] = pow(10,19) ;// initial n
floatinphys[30] = 0.0 ;// initial v 
floatinphys[31] = 100.0 ;// initial T 
floatinphys[32] = 6.2*pow(10,17) ;// initial a
floatinphys[33] = 0.0 ;// initial vn
floatinphys[34] = 5.5*pow(10,17) ;// initial m
floatinphys[35] = 2.7*pow(10,18) ;// initial ncor 
floatinphys[36] = 1.0*pow(10,14) ;// initial core neutral

floatinphys[39] = 6.28*pow(10,17);
floatinphys[40] = 6.28*pow(10,17);
floatinphys[41] = 5.4*pow(10,16); 
floatinphys[42] = 6.28*pow(10,17);
floatinphys[43] = 6.28*pow(10,17); // initial nb(5)

floatinphys[44] = 3.89*pow(10,18);
floatinphys[45] = 6.44*pow(10,18);
floatinphys[46] = 1.87*pow(10,17);
floatinphys[47] = 6.44*pow(10,18);
floatinphys[48] = 3.89*pow(10,18);// initial mb(5)

// limits
floatinphys[49] = pow(10,10);// min den
floatinphys[50] = pow(10,25);// max den
floatinphys[51] = 0.1000;// min tem

// sol behavior
floatinphys[59] = 6 ;// gamma
floatinphys[60] = 3.3436*pow(10,-27) ;// mass
floatinphys[61] = 30 ;// energy loss ion
floatinphys[62] = 0.4 ;// recycling
floatinphys[63] = 1.0;// mol_rec
floatinphys[64] = 5.0*pow(10,-4) ;// neutral residence time
floatinphys[65] = 8.0*pow(10,-3) ;// molecule residence time
floatinphys[66] = 0.0000 ; // far sol ion losses 
// new array
floatinphys[69] =2.0*pow(10,-3); 
floatinphys[70] = 0.0*pow(10,-12);
floatinphys[71] = 0.0*pow(10,-12);
floatinphys[72] = 0.0*pow(10,-12);
floatinphys[73] = 0.0*pow(10,-12);  // imp con 1:5

// reservoirs
floatinphys[79] = pow(10,-22);
floatinphys[80] = pow(10,-22);
floatinphys[81] = 2.05*pow(10,2);
floatinphys[82] = pow(10,-22); 
floatinphys[83] = pow(10,-22);// cor res ato pmp (1:5)

floatinphys[84] = pow(10,-22);
floatinphys[85] = pow(10,-22); 
floatinphys[86] = 2.05*pow(10,2); 
floatinphys[87] = pow(10,-22);
floatinphys[88] = pow(10,-22); // cor res mol pmp(1:5)

floatinphys[89] = 1.0599; 
floatinphys[90] = 1.826200;
floatinphys[91] = 30.000;
floatinphys[92] = 1.8268;
floatinphys[93] = 1.0599;// ext vol 1:5

floatinphys[99] = 0.000;
floatinphys[100] = 0.000;
floatinphys[101] = 6.67*pow(10,2);// ext ato ext 1:3
floatinphys[102] = 0.000;
floatinphys[103] = 0.000; 
floatinphys[104] = 4.717*pow(10,2) ;// ext mol ext 1:3

floatinphys[109] = 0.000;
floatinphys[110] = 0.000;
floatinphys[111] = 1.25*pow(10,4);
floatinphys[112] = 6.25*pow(10,2);
floatinphys[113] = 1.875*pow(10,3);// pmp rat n(1:5)

floatinphys[114] = 0.000;
floatinphys[115] = 0.000; 
floatinphys[116] = 1.25*pow(10,4);
floatinphys[117] = 6.25*pow(10,2); 
floatinphys[118] = 1.875*pow(10,3);// pmp rat m(1:5)

// core
floatinphys[129] = 0.2000 ;// core conf tim
floatinphys[130] = 0.0*pow(10,-22);// cor sol ato exchange
floatinphys[131] = 0.0*pow(10,-22);// cor sol mol exchange
floatinphys[132] = 4.877; // cor volumes


 int call_from_extern = 1;
    // store adress of settings in pointers
int *call_from_externPtr = &call_from_extern; 
//print("C size floatinphys %4.1d",size(floatinphys));
 // grid data
 #ifdef INITIALIZE_DIV1D
 int Nx = intinnum[0];
 int *NxPtr = &Nx;
 int init_grid_fortran = 1;
 int init_prof_fortran = 1;
 int *init_grid_fortranPtr = &init_grid_fortran;   
 int *init_prof_fortranPtr = &init_prof_fortran;
	
// arrays are automatically pointers, but values are not :)

 // grid data 
 double E_x[Nx]; 
 double E_xcb[Nx+1];
 double E_delta_x[Nx];
 double E_delta_xcb[Nx+1];
 double E_B_field[Nx];
 double E_B_field_cb[Nx+1];
 double E_B_trans[Nx];
 double E_B_trans_cb[Nx+1];
 double E_R_cc[Nx], E_R_cb[Nx+1], E_Area_extern[Nx], E_sintheta_cc[Nx], E_sintheta_cb[Nx+1];
 double E_sol_width_pol[Nx], E_sol_width_pol_cb[Nx+1], E_volumes[Nx];
 double E_gas_puff_profile[Nx], E_core_source_profile_q[Nx], E_core_source_profile_n[Nx];
 int E_i_omp, E_i_Xpoint[2]; 
 int *E_i_ompPtr = &E_i_omp;
 // int *E_i_XpointPtr = &E_i_Xpoint; // arrays are already pointer
 double E_mid_point, E_X_omp[2]; 
 double *E_mid_pointPtr = &E_mid_point;
 int E_i_baffle[2];
 //double 
 // mid_point also a pointer
			  // E_R_cc, E_R_cb, E_Area_extern, E_sintheta_cc, E_sintheta_cb, E_sol_width_pol, E_sol_width_pol_cb, E_volumes,
			  // E_gas_puff_profile, E_core_source_profile,

 // initial plasma condittions
 double E_density[Nx];
 double E_velocity[Nx];
 double E_temperature[Nx];
 double E_neutral[Nx];
 double E_neutral_velocity[Nx];
 double E_molecule[Nx];
 double E_Gamma_n[Nx+1]; 
 double E_Gamma_mom[Nx+1];
 double E_q_parallel[Nx+1];
 double E_Gamma_neutral[Nx+1];
 double E_Gamma_mom_neutral[Nx+1];
 double E_Gamma_molecule[Nx+1];
 double E_Source_n[Nx];
 double E_Source_v[Nx];
 double E_Source_Q[Nx];
 double E_Source_neutral[Nx];
 double E_Source_vn[Nx];
 double E_Source_molecule[Nx];

 double E_core2sol_flux[Nx];
 double E_extern2sol_flux[Nx];
 double E_core2sol_mol[Nx];
 double E_extern2sol_mol[Nx];
 double E_sol2extern_ion_flux[Nx];

 for(int i=0; i < Nx; ++i )
 {
 E_density[i] = pow(10,20);
 E_velocity[i] = 0.0;
 E_temperature[i]= 20; 
 E_neutral[i] = pow(10,18);
 E_neutral_velocity[i]= 0.0;
 E_molecule[i] = pow(10,16);
 E_Source_n[i]=0.0;
 E_Source_v[i]=0.0;
 E_Source_Q[i]=0.0;
 E_Source_neutral[i]=0.0;
 E_Source_vn[i]=0.0;
 E_Source_molecule[i]=0.0;
 E_extern2sol_flux[i] = 0.0;  
 E_extern2sol_mol[i] = 0.0;
 E_core2sol_flux[i] = 0.0;
 E_core2sol_mol[i] = 0.0;
 E_sol2extern_ion_flux[i]=0.0;
 } 
for(int i=0; i< Nx+1; ++i)
{
 E_Gamma_n[i] = 0.0; 
 E_Gamma_mom[i] =0.0;
 E_q_parallel[i]=0.0;
 E_Gamma_neutral[i]=0.0;
 E_Gamma_mom_neutral[i]=0.0;
 E_Gamma_molecule[i] = 0.0;
}


double E_extern_neutral_flux[3];
double E_extern_molecule_flux[3];
double E_sol2extern_flux[5];
double E_sol2extern_mol[5];
double E_tar2extern_flux[5]; 		
double E_tar2extern_mol[5]; 		     
double E_Source_extern[5];
double E_Source_extern_mol[5];
double E_extern2core_flux[5];
double E_extern2core_mol[5];

/*for(int i=0; i< 5; ++i)
{
E_Source_extern[i] = 0.0;
E_Source_extern_mol[i] = 0.0;
}*/

double zero_init1 = 0.00;
double zero_init2 = 0.00;
double zero_init3 = 0.00;
double zero_init4 = 0.00;

double *E_sum_sol2extern_ion_flux;
E_sum_sol2extern_ion_flux = &zero_init1;
double *E_sum_sol2extern_ion_mol;
E_sum_sol2extern_ion_mol = &zero_init2;
double *E_sol2core_flux;
E_sol2core_flux = &zero_init3;
double *E_sol2core_mol;
E_sol2core_mol = &zero_init4;

double E_Gamma_core2sol_val = pow(10,21);
double *E_Gamma_core2sol; // initial value? 10^20
E_Gamma_core2sol = &E_Gamma_core2sol_val;
double E_Source_core_val = 0.0; 
double *E_Source_core;
E_Source_core = &E_Source_core_val;
double E_molecule_pump[5];
double E_neutral_pump[5];

double E_core_density_val = 2.7*pow(10,18); // floatinphys ... 
double *E_core_density;
E_core_density = &E_core_density_val; // initial value // referencing pointer to E_core_density_val, 
//E_core_density_val = *E_core_density; // dereferencing pointer value to core-density-val
double E_extern_neutral_density[5]; // initial value
double E_extern_molecule_density[5]; // initial value


E_extern_neutral_density[0] = 6.28*pow(10,17);
E_extern_neutral_density[1] = 6.28*pow(10,17);
E_extern_neutral_density[2] = 5.4*pow(10,16); 
E_extern_neutral_density[3] = 6.28*pow(10,17);
E_extern_neutral_density[4] = 6.28*pow(10,17); // initial nb(5)

E_extern_molecule_density[0] = 3.89*pow(10,18);
E_extern_molecule_density[1] = 6.44*pow(10,18);
E_extern_molecule_density[2] = 1.87*pow(10,17);
E_extern_molecule_density[3] = 6.44*pow(10,18);
E_extern_molecule_density[4] = 3.89*pow(10,18);// initial mb(5)


  //#ifdef RUN_DIV1D  
 // running parameters
 double start_time = 0.0;
 double *start_timePtr = &start_time;

 double end_time = 0.0;
 double *end_timePtr = &end_time;

 int nout = intinnum[1];
 int *noutPtr = &nout; 

 double delta_t = floatinnum[60];
 double *delta_tPtr = &delta_t; // this pointer holds the adress of delta_t 


 double E_imp_con[nout][5];
 //double *E_imp_conPtr = &E_imp_con; // pointer types for matrix might be incompatible
 double E_neu[nout];
 double E_dneu[nout];
 double E_nb[nout][5];
 double E_mb[nout][5];
 double E_gas[nout];
 double E_rec[nout];
 double E_qpar_x[nout];
 double E_red_frc[nout];
 double E_Q_core[nout];
 double E_Gamma_core[nout];
 double E_core_neutral_density[nout];
 double E_neutral_puff[nout][5];
 double E_molecule_puff[nout][5];
 double E_core_fuelling[nout];

 for(int i=0; i < nout; ++i )
 {
 // this seems to be the right orientation.
 // C sends imp_con[nout][4]
 // F reads imp_con[4][nout]
 E_imp_con[i][0] = 0.03; //,0.0,0.0,0.0,0.0};
 E_imp_con[i][1] = 0.0;
 E_imp_con[i][2] = 0.0;
 E_imp_con[i][3] = 0.0;
 E_imp_con[i][4] = 0.0;
 E_neu[i] = pow(10,20);
 E_dneu[i] = 0.0;
 E_nb[i][0] = pow(10,18);
 E_nb[i][1] = pow(10,18);
 E_nb[i][2] = pow(10,18);
 E_nb[i][3] = pow(10,18);
 E_nb[i][4] = pow(10,18);
 E_mb[i][0] = pow(10,16);
 E_mb[i][1] = pow(10,16);
 E_mb[i][2] = pow(10,16);
 E_mb[i][3] = pow(10,16);
 E_mb[i][4] = pow(10,16);
 E_gas[i] = 0.0;
 E_rec[i] = 0.99;
 E_qpar_x[i] = 20*pow(10,6);
 E_red_frc[i] = 0.0;
 E_Q_core[i] = 500.0*pow(10,3);
 E_Gamma_core[i] = 1.0*pow(10,20);
 E_core_neutral_density[i]= pow(10,14); 
/*PUFF_RATE_MOLECULE =  0.00000000000000E+00  0.00000000000000E+00  4.34277978975998E+21
  6.31896200000000E+21  4.10031500000000E+21*/
/* PUFF_RATE_NEUTRAL =  0.00000000000000E+00  0.00000000000000E+00  7.11753051696510E+20
  9.40378875000000E+19  3.82575800000000E+20 */
 E_neutral_puff[i][0] = 0.0*pow(10,-20);
E_neutral_puff[i][1] = 0.0*pow(10,-20);
E_neutral_puff[i][2] = 7.1*pow(10,20);
E_neutral_puff[i][3] = 9.4*pow(10,19);
E_neutral_puff[i][4] = 3.8*pow(10,20);
  E_molecule_puff[i][0] = 0.0*pow(10,-20);
  E_molecule_puff[i][1] = 0.0*pow(10,-20);
  E_molecule_puff[i][2] = 4.3*pow(10,21);
  E_molecule_puff[i][3] = 6.3*pow(10,21);
  E_molecule_puff[i][4] = 4.1*pow(10,21);
E_core_fuelling[i] = 0.0;
 }

 printf("C number of cells Nx =  %d \n", Nx); 
 printf("C Initial conditions: \n");
 printf("C Density at i=3 = %4.1e \n", E_density[4]);
 printf("C velocity at i = 7 = %4.1e \n", E_velocity[8]);
 printf("C temperature at i = 29 = %4.1e \n", E_temperature[30]);
 printf("C neutral at i = 299 = %4.1e \n",E_neutral[300]);
 for(int i=0; i<nout; ++i) {
    printf("C impurities %4.1e \n",E_imp_con[0][i]);
 }

 //#endif
 #endif
 #endif
    
 #ifdef INITIALIZE_SETTINGS
 printf("C Initialize DIV1D settings \n");
 initialize_div1d_settings_(floatinnum, intinnum, loginnum, floatinphys, intinphys, call_from_externPtr);
 #ifdef INITIALIZE_DIV1D
 printf("C Initializing External Library \n");
 printf("C pushed information in xgrid at i = 3 = %4.1e \n",E_x[4]);
 printf("C pushed information in density at i = 199 = %4.1e \n",E_density[200]);   
 printf("C pushed information in B_field at i = 199 = %4.1e \n",E_B_trans[200]);   

 initialize_div1d_arrays_(E_density, E_velocity, E_temperature, E_neutral, E_neutral_velocity, E_molecule,
		           E_x, E_xcb, E_delta_x, E_delta_xcb, E_B_field, E_B_field_cb, E_B_trans, E_B_trans_cb,
			   E_R_cc, E_R_cb, E_Area_extern, E_sintheta_cc, E_sintheta_cb, E_sol_width_pol, E_sol_width_pol_cb, E_volumes,
			   E_gas_puff_profile, E_core_source_profile_q, E_core_source_profile_n,
		           init_grid_fortranPtr, init_prof_fortranPtr, NxPtr, E_i_ompPtr, E_i_Xpoint, E_i_baffle, E_mid_pointPtr, E_X_omp);

 printf("C received information in xgrid at i = 3 = %4.1e \n",E_x[4]);
 printf("C received information in density at i = 199 = %4.1e \n",E_density[200]);
 printf("C Initialized DIV1D \n");

 #ifdef RUN_DIV1D
printf("C Calling External Library \n ");
/* run_div1d_(E_density, E_velocity, E_temperature, E_neutral, 
            E_Gamma_n, E_Gamma_mom, E_q_parallel, E_neutral_flux, // Output
            E_extern_neutral_flux, E_sol2extern_flux, E_extern2sol_flux, E_tar2extern_flux, // Output
            E_sol2core_flux, E_core2sol_flux,  // Output
            E_Source_n, E_Source_v, E_Source_Q, E_Source_neutral, E_Source_extern, // Output
            start_timePtr, end_timePtr, noutPtr, delta_tPtr, 
            E_imp_con, E_neu, E_dneu, E_ngb, E_gas, E_rec,  E_qpar_x, E_red_frc, E_Q_core, E_Gamma_core); */
 printf("C Core vars %4.1e %4.1e %4.1e %4.1e \n %4.1e ", E_core_density, E_Gamma_core2sol,E_sol2core_flux, E_sol2core_mol, E_Source_core);

// step a few timesteps
for(int i=0; i<10; ++i) {
printf("C start step %d, at time %4.1e \n" , i, start_time);

run_div1d_(E_density,E_velocity, E_temperature, E_neutral, E_neutral_velocity, E_molecule, //6
          E_Gamma_n, E_Gamma_mom, E_q_parallel, E_Gamma_neutral, E_Gamma_mom_neutral, E_Gamma_molecule, //12
	  E_core2sol_flux, E_core2sol_mol, E_sol2extern_ion_flux,  //15
	  E_Source_n, E_Source_v, E_Source_Q, E_Source_neutral, E_Source_vn, E_Source_molecule, //21
	  E_extern_neutral_density, //22
	  E_extern_neutral_flux, E_sol2extern_flux, E_extern2sol_flux, E_tar2extern_flux, //26
	  E_extern2core_flux, E_sum_sol2extern_ion_flux, //28
	  E_Source_extern, E_neutral_pump, //30
	  E_extern_molecule_density,  //31
	  E_extern_molecule_flux, E_sol2extern_mol, E_extern2sol_mol, E_tar2extern_mol, //35
	  E_extern2core_mol, E_sum_sol2extern_ion_mol,  //37
	  E_Source_extern_mol, E_molecule_pump, //39
	  E_core_density, E_Gamma_core2sol,E_sol2core_flux, E_sol2core_mol, E_Source_core,//44
          start_timePtr, end_timePtr, noutPtr, delta_tPtr, //48
	  E_imp_con, E_neu, E_dneu, E_nb, E_mb, E_gas, E_rec, E_qpar_x, E_red_frc, //57
          E_Q_core, E_Gamma_core, E_core_neutral_density, // 60
	  E_neutral_puff, E_molecule_puff, E_core_fuelling); //63

start_time = *end_timePtr;
printf("C finished %d, at time %4.1e \n" , i, *end_timePtr);
start_timePtr = &start_time;
}

/* run_div1d_(E_density, E_velocity, E_temperature, E_neutral,
		 E_Gamma_n, E_Gamma_mom, E_q_parallel, E_neutral_flux,
		 E_Source_n, E_Source_v, E_Source_Q, E_source_neutral,
		 start_timePtr, end_timePtr, noutPtr, delta_tPtr,
		 E_imp_con, E_neu, E_dneu, E_ngb, E_gas, E_rec, E_qpar_x, E_red_frc);  // , E_Q_core, E_Gamma_core */

 printf("C received information in density at i = 199 = %4.1e \n",E_density[200]);
 printf("C received information in extern2sol_flux at i = 199 = %4.1e \n",E_extern2sol_flux[200]);
 printf("C finished nout steps with DIV1D\n");
 #endif
 #endif
 #endif
 printf("C Done testing the DIV1D testscript from C\n");

 return 0;

}
