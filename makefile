# makefile for div1d programme

# Fortran compiler used
FC = ifort         # for INTEL fortran compiler
# FC = gfortran      # for GNU   fortran compiler

# compiler options
FOPT_DVODE =  -module obj -g -O3 # for intel fortran compiler
# FOPT_DVODE = -p -module obj -g -O3 # for intel fortran compiler
# FOPT_DVODE =  -module obj -g -check all -debug all # for intel fortran compiler with all checks and debugging info
# FOPT_DVODE =  -Jobj -ffree-line-length-none -g -fbacktrace # for GNU fortran compiler with all checks and debugging info
FOPT =  -module obj -g -O3 # for intel fortran compiler
# FOPT =  -module obj -g -check all -debug all # for intel fortran compiler with all checks and debugging info
# FOPT =  -Jobj -ffree-line-length-none -g -fbacktrace # for GNU fortran compiler with all checks and debugging info
# FOPT = -p -module obj -g -O3 # for intel fortran compiler and profiling with gprof

OBJECTS = obj/constants.o \
          obj/physics_routines.o \
          obj/grid_data.o \
          obj/numerics_parameters.o \
          obj/physics_parameters.o \
          obj/plasma_data.o \
          obj/reaction_rates.o \
          obj/interpolation.o \
          obj/dvode_f90_m.o \
          obj/dlsode.o \
          obj/opkda1.o \
          obj/opkda2.o \
          obj/rk4.o


obj/div1d.exe : src/div1d_main.f90 $(OBJECTS)
	$(FC) src/div1d_main.f90 $(OBJECTS)  $(FOPT) -o obj/div1d.exe


obj/constants.o : src/constants.f90
	$(FC) src/constants.f90 $(FOPT) -c -o obj/constants.o


obj/grid_data.o : src/grid_data.f90\
                   obj/numerics_parameters.o \
                   obj/physics_parameters.o
	$(FC) src/grid_data.f90 $(FOPT) -c -o obj/grid_data.o


obj/numerics_parameters.o : src/numerics_parameters.f90
	$(FC) src/numerics_parameters.f90 $(FOPT) -c -o obj/numerics_parameters.o


obj/physics_parameters.o : src/physics_parameters.f90
	$(FC) src/physics_parameters.f90 $(FOPT) -c -o obj/physics_parameters.o


obj/physics_routines.o : src/physics_routines.f90\
                   obj/constants.o \
                   obj/grid_data.o \
                   obj/reaction_rates.o \
                   obj/physics_parameters.o
	$(FC) src/physics_routines.f90 $(FOPT) -c -o obj/physics_routines.o


obj/plasma_data.o : src/plasma_data.f90\
                   obj/numerics_parameters.o \
                   obj/physics_parameters.o \
                   obj/physics_routines.o \
                   obj/interpolation.o
	$(FC) src/plasma_data.f90 $(FOPT) -c -o obj/plasma_data.o


obj/reaction_rates.o : src/reaction_rates.f90
	$(FC) src/reaction_rates.f90 $(FOPT) -c -o obj/reaction_rates.o


obj/interpolation.o : src/interpolation.f90
	$(FC) src/interpolation.f90 $(FOPT) -c -o obj/interpolation.o


obj/dvode_f90_m.o : src/dvode_f90_m.f90
	$(FC) src/dvode_f90_m.f90 $(FOPT_DVODE) -c -o obj/dvode_f90_m.o


obj/dlsode.o : src/dlsode.f\
                   obj/opkda1.o \
                   obj/opkda2.o
	$(FC) src/dlsode.f $(FOPT) -c -o obj/dlsode.o


obj/opkda1.o : src/opkda1.f
	$(FC) src/opkda1.f $(FOPT) -c -o obj/opkda1.o


obj/opkda2.o : src/opkda2.f
	$(FC) src/opkda2.f $(FOPT) -c -o obj/opkda2.o


obj/rk4.o : src/rk4.f90
	$(FC) src/rk4.f90 $(FOPT) -c -o obj/rk4.o



clean:
	rm obj/*.o obj/*.mod obj/div1d.exe
