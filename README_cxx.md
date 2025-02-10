Shallow water model
===================

The shallow water model is one of the simplest examples of an atmosphere or ocean model. 

Methods to optimise
* Compiler flags
* Fortran vs C++ 
* 2D array vs 1D array 
* OpenMP parallelisation `OMP_NUM_THREADS=4` 
* Pre-compute constants e.g. `1/dx` etc
* Write time-averages rather than each step
* Use vectorisation `-march=native` 
* Single precision vs double precision
* Turn off non-linear terms
