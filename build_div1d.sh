#!/bin/bash
PWD=$(pwd)

 make clean
# or
rm -r obj
mkdir obj

make

chmod 774 $PWD"/obj/div1d.exe"

