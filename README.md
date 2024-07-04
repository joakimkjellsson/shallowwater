1-layer linear shallow water model
==================================

The model integrates the linearized shallow water equations (ref. equation 3.112 in "Atmospheric and Oceanic Fluid Dynamics" by Vallis, 4th Ed.) on a C-grid. The first time step is Euler forward (using information in step 1 to get to step 2). The following steps are done by a leap frog scheme (using information in step n-1 and n to get to n+1).
The "nt" time steps are divided into "nsubcycles" steps. Thus, the model computes nt*nsubcycles steps, but only saves nt steps.

This model is intended as a minimal example of an atmosphere or ocean model. It is written in Fortran 90 with netCDF support. 

Prerequisites
-------------

Install a Fortran compiler (Mac OS): 
1) Download MacPorts and follow instructions to install. 
2) Then install GCC. In this case we use GCC v13 (but any later version should be fine too): 
```bash
sudo port install gcc13
```
If you chose to install MacPorts in `/opt/local/` you should now have `/opt/local/bin/gfortran-mp-13` (assuming you installed v13). 

3) Install netcdf (optional): 
```bash
sudo port install netcdf netcdf-fortran
```
You should now have the netCDF Fortran libraries in `/opt/local/lib/`. You should set `NCDF_ROOT=/opt/local/` in the `makefile.netcdf`. 

Install a Fortran compiler (Linux): 
1a) Use `apt-get` (Ubuntu, Linux Mint, etc)
```bash
sudo apt-get install gfortran
```
1b) Use `yum` (Fedora, Red Hat etc)
```bash
sudo yum install gcc-gfortran
```
2) Install netcdf (optional)
```bash
sudo apt-get install netcdf netcdf-fortran
``` 
(or similar with `yum`)

Install Fortran compiler (Windows):
1) Install Linux and follow instructions above ;-) 

Compile the model
-----------------

If you want netCDF support use the `makefile.netcdf`, otherwise use `makefile.binout`. 
Which ever you choose you will need to set the following: 
```bash
FC = gfortran (or whatever your compiler name is)
NCDF_ROOT = /opt/local/ (or whereever your netcdf is installed)
```
Note: You can usually find the netcdf directory with the command
```bash
nf-config --prefix
```

Now compile with `make -f makefile.netcdf` (for netcdf) or `make -f makefile.binout` (without netcdf). 

Running the model
-----------------

To run the lab3 experiment (geostrophically balanced flow on a beta plane): 

Link the namelist to a file `namelist`
```bash
ln namelist_lab3 namelist
```

Run the model
```bash
./shallowwater_lab3.x
```

Changing model parameters
-------------------------

Model parameters are controlled via namelists, which are text files that Fortran reads. 
The model will always read the file `namelist`, which can be a link to another file, e.g. `namelist_lab1`. 
The model sets a lot of variables by default, e.g. domain size, Coriolis parameter etc, but the namelist overwrites the model defaults. 
The model parameters should therefore be changed by changing the namelist, not the model code. 

Output
------

Results will be saved to `u_test.nc`, `v_test.nc` and `h_test.nc` netCDF files. 
