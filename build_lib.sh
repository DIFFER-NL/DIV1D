#!/usr/bin/env bash
# makes the div1d library shared object file

PWD=$(pwd)


cd obj

# compile dynamic library
ifort -shared -fPIC -v   -o  libdiv1d.so *.o 

nm -D libdiv1d.so # checks the parameters and functions that can be accessed 
# for information on the nm command look at linux nm command @thegeekstuff.com

cd ..



