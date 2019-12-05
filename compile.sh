#!/bin/bash

fortran_input="./src/main/fortran/square_cube.f90"
fortran_obj="./squarecube.o"
gfortran-9 -fPIC -o $fortran_obj -c $fortran_input
gcc -shared -o libfoo.so squarecube.o

c_src="./src/jni/c/square_cube.c"
c_obj="./test"

gcc -L/home/jamie/sqslp_in_java -L/usr/lib/gcc/x86_64-linux-gnu/9 -o $c_obj $c_src -lfoo -lgfortran
