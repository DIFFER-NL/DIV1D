# makefile for div1d programme

# Fortran compiler used
FC = ifort         # for INTEL fortran compiler

# compiler options
FOPT =  -module obj -g  # for intel fortran compiler

OBJECTS = obj/constants.o \
          obj/physics_routines.o \
          obj/grid_data.o \
          obj/numerics_parameters.o \
          obj/physics_parameters.o \
          obj/plasma_data.o \
          obj/dvode_f90_m.o \
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
                   obj/physics_parameters.o
	$(FC) src/physics_routines.f90 $(FOPT) -c -o obj/physics_routines.o


obj/plasma_data.o : src/plasma_data.f90\
                   obj/numerics_parameters.o \
                   obj/physics_parameters.o \
                   obj/physics_routines.o
	$(FC) src/plasma_data.f90 $(FOPT) -c -o obj/plasma_data.o


obj/dvode_f90_m.o : src/dvode_f90_m.f90
	$(FC) src/dvode_f90_m.f90 $(FOPT) -c -o obj/dvode_f90_m.o


obj/rk4.o : src/rk4.f90
	$(FC) src/rk4.f90 $(FOPT) -c -o obj/rk4.o



clean:
	rm obj/*.o obj/*.mod obj/div1d.exe
