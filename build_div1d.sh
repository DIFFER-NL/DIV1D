#!/bin/bash
PWD=$(pwd)

mkdir obj

make

chmod 774 $PWD"/obj/div1d.exe"

