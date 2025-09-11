# makefile for div1d programme

# Fortran compiler used
FC = ifort         # for INTEL fortran compiler
# FC = gfortran      # for GNU   fortran compiler

# compiler options

 FOPT_DVODE = -p -module obj -g -O3 -fPIC # for intel fortran compiler
# FOPT_DVODE =  -module obj -g -fPIC -check all -debug all # for intel fortran compiler with all checks and debugging info
# FOPT_DVODE =  -Jobj -ffree-line-length-none -g -fbacktrace # for GNU fortran compiler with all checks and debugging info

 FOPT =  -module obj -g -fPIC -O3 # for intel fortran compiler
# FOPT =  -module obj -g -fPIC -check all -debug all # for intel fortran compiler with all checks and debugging info
# FOPT =  -Jobj -ffree-line-length-none -g -fbacktrace # for GNU fortran compiler with all checks and debugging info
# FOPT = -p -module obj -g -O3 # for intel fortran compiler and profiling with gprof

#LOPT = -module obj -g -03  -fPIC
LOPT = -module obj -g fPIC -check all -debug all

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
          obj/rk4.o\
	  obj/experiments.o\
	  obj/radiative_cooling_functions_post1977.o\
	  obj/div1d_step.o

LIBOBJ = obj/constants.o\
	 obj/physics_routines.o\
	 obj/grid_data.o\
	 obj/numerics_parameters.o\
	 obj/physics_parameters.o\
	 obj/plasma_data.o \
         obj/reaction_rates.o \
	 obj/interpolation.o \
         obj/dvode_f90_m.o \
	 obj/dlsode.o \
         obj/opkda1.o \
	 obj/opkda2.o \
	 obj/rk4.o\
	 obj/experiments.o\
	 obj/radiative_cooling_functions_post1977.o   

obj/div1d.exe : src/div1d_main.f90 $(OBJECTS)
	$(FC) src/div1d_main.f90 $(OBJECTS)  $(FOPT) -o obj/div1d.exe

# dynamic library
#obj/div1dlib.so : src/div1d_solve.f90 $(LIBOBJ)
#	$(FC)  -o div1dlib.so -G -fPIC -ztext -hdiv1dlib.so ./../src/div1d_solve.f90

#obj/div1dlib.so : src/div1d_solve.f90 $(LIBOBJ)
#	$(FC) src/div1d_solve.90 $(LIBOBJ) $(LOPT) -o obj/div1dlib.so

# static library
#ar cr obj/div1dlib.a obj/div1d_solve.o

#obj/div1d_test.exe : src/div1d_test.f90 $(OBJECTS)
#	$(FC) src/div1d_test.f90 $(OBJECTS)  $(FOPT) -o obj/div1d_test.exe

#obj/div1d_solve.so : src/div1d_solve.f90 $(LIBOBJ)    
#	$(FC) src/div1d_solve.f90 -I$(LIBOBJ)  $(LOPT)  -o obj/div1d_solve.so -G  -ztext -hdiv1dsolve.so 

obj/div1d_step.o : src/div1d_step.f90 $(LIBOBJ)    		    
	$(FC) src/div1d_step.f90  $(FOPT) -c -o obj/div1d_step.o

#obj/div1d_step.o : src/div1d_step.f90\
#		  obj/constants.o\
#		  obj/physics_routines.o\
#		  obj/grid_data.o\
#		  obj/numerics_parameters.o\
#		  obj/physics_parameters.o\
#	          obj/plasma_data.o \
 #       	  obj/reaction_rates.o \
#	          obj/interpolation.o \
 #         	  obj/dvode_f90_m.o \
#	          obj/dlsode.o \
 #       	  obj/opkda1.o \
#	          obj/opkda2.o \
#	          obj/rk4.o\
#		  obj/experiments.o\
#		  obj/radiative_cooling_functions_post1977.o      		    
#	$(FC) src/div1d_step.f90 $(FOPT) -c -o obj/div1d_step.o


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


obj/dvode_f90_m.o : src/dvode_f90_m.f90
	$(FC) src/dvode_f90_m.f90 $(FOPT_DVODE) -c -o obj/dvode_f90_m.o

obj/physics_routines.o : src/physics_routines.f90\
                   obj/constants.o \
                   obj/grid_data.o \
                   obj/reaction_rates.o \
                   obj/physics_parameters.o \
                   obj/experiments.o \
		   obj/dvode_f90_m.o
	$(FC) src/physics_routines.f90 $(FOPT) -c -o obj/physics_routines.o


obj/plasma_data.o : src/plasma_data.f90\
                   obj/numerics_parameters.o \
                   obj/physics_parameters.o \
                   obj/physics_routines.o \
                   obj/interpolation.o
	$(FC) src/plasma_data.f90 $(FOPT) -c -o obj/plasma_data.o


obj/reaction_rates.o : src/reaction_rates.f90\
		   obj/radiative_cooling_functions_post1977.o
	$(FC) src/reaction_rates.f90 $(FOPT) -c -o obj/reaction_rates.o


obj/radiative_cooling_functions_post1977.o : src/radiative_cooling_functions_post1977.f90
	$(FC) src/radiative_cooling_functions_post1977.f90 $(FOPT) -c -o obj/radiative_cooling_functions_post1977.o


obj/interpolation.o : src/interpolation.f90
	$(FC) src/interpolation.f90 $(FOPT) -c -o obj/interpolation.o



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


obj/experiments.o : src/experiments.f90
	$(FC) src/experiments.f90 $(FOPT) -c -o obj/experiments.o

clean:
	rm obj/*.o obj/*.mod obj/div1d.exe obj/div1d_test.exe

#Commands below are related to test automation. 
TEST_BUILD_DIR=build
TEST_TYPE ?=
TEST_ARGS ?= --v

.PHONY: test

test:
	@cd test && \
	if [ ! -d "$(TEST_BUILD_DIR)" ]; then mkdir -p "$(TEST_BUILD_DIR)"; fi && \
	cd "$(TEST_BUILD_DIR)" && FC=$(FC) cmake .. && make
	@FAILED=0; \
	cd test/"$(TEST_BUILD_DIR)"; \
	if [ -z "$(TEST_TYPE)" ]; then \
		if [ -d "unit_tests" ]; then \
			for f in unit_tests/run_*_tests; do \
				if [ -x "$$f" ]; then \
					echo ""; \
					echo "\033[1;33m======= Running $$f =======\033[0m"; \
					echo ""; \
					cp "$$f" . && \
					./$$(basename $$f) $(TEST_ARGS); \
					STATUS=$$?; \
					rm -f ./$$(basename $$f); \
					if [ $$STATUS -ne 0 ]; then \
						echo "\033[1;31m $$f failed \033[0m"; \
						FAILED=1; \
					fi; \
				fi; \
			done; \
		fi; \
		for dir in *_tests; do \
			if [ "$$dir" = "unit_tests" ]; then continue; fi; \
			if [ -d "$$dir" ]; then \
				for f in "$$dir"/run_*_tests; do \
					if [ -x "$$f" ]; then \
						echo ""; \
						echo "\033[1;33m======= Running $$f =======\033[0m"; \
						echo ""; \
						cp "$$f" . && \
						./$$(basename $$f) $(TEST_ARGS); \
						STATUS=$$?; \
						rm -f ./$$(basename $$f); \
						if [ $$STATUS -ne 0 ]; then \
							echo "\033[1;31m $$f failed \033[0m"; \
							FAILED=1; \
						fi; \
					fi; \
				done; \
			fi; \
		done; \
	else \
		for dir in $(TEST_TYPE)_tests; do \
			if [ -d "$$dir" ]; then \
				for f in "$$dir"/run_*_tests; do \
					if [ -x "$$f" ]; then \
						echo ""; \
						echo "\033[1;33m======= Running $$f =======\033[0m"; \
						echo ""; \
						cp "$$f" . && \
						./$$(basename $$f) $(TEST_ARGS); \
						STATUS=$$?; \
						rm -f ./$$(basename $$f); \
						if [ $$STATUS -ne 0 ]; then \
							echo "\033[1;31m $$f failed \033[0m"; \
							FAILED=1; \
						fi; \
					fi; \
				done; \
			fi; \
		done; \
	fi; \
    mkdir -p test_results && \
    find . -maxdepth 1 -name '*.xml' -exec mv {} test_results/ \; && \
	exit $$FAILED

