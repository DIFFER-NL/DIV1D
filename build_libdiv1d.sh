#!/usr/bin/env bash
# makes the div1d library

PWD=$(pwd)
./build_div1d.sh

cd obj

# compile dynamic library
ifort -shared -fPIC -v   -o  libdiv1d.so *.o 

nm -D libdiv1d.so # checks the parameters and functions that can be accessed 
# for information on the nm command look at linux nm command @thegeekstuff.com

cd ..

# linking the C test wrapper to the shared library
gcc -L$PWD/obj -o ./obj/test_libdiv1d.exe $PWD/src/test_libdiv1d.c -I$PWD/src -ldiv1d -ldl -lstdc++

# export library path to find it during execution
export LD_LIBRARY_PATH=$PWD/obj:$LD_LIBRARY_PATH

# call the library
./obj/test_libdiv1d.exe



# legacy
# compile static library (presently not used and not tested)
# ar cr libdiv1d.a div1d_solve.o
# compile the C function that calls the library
# gcc text_libdiv1d.c div1d_solve.o -g -v
# gcc -v -ldl test_libdiv1d.c -o $PWD/obj/test_libdiv1d_c.exe

