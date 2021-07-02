#!/bin/bash

rm *.mod

# optimizations: -O3 -ipo (very important!!)
export F_UFMTENDIAN=big
mpif90 -check all mpi.f90 global.f90 mass_flow.f90 interpolation.f90 equations.f90 boundary_conditions.f90 pressure.f90 input_output.f90 statistics.f90 initialization.f90 eqwm.f90 subgrid.f90 wallmodel.f90 finalization.f90 projection.f90 time_integration.f90 monitor.f90  main.f90 -o channel_WM_large_0.25 -lfftw3_mpi -lfftw3 /usr/lib64/liblapack.so.3.2.1 -I/share/apps/lib/fftw/3.3.4/mvapich2-2.0rc1-intel-14/include/ -L/share/apps/lib/fftw/3.3.4/mvapich2-2.0rc1-intel-14/lib/

export F_UFMTENDIAN=big

