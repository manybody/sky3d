Sky3d is a nuclear time-dependent Hartree-Fock code.
It represents wave functions in a symmetry-unrestricted three-dimensional cartesian box, and operates in static and time-dependent modes.  The static mode is mainly for the generation of initial conditions for time-depdendent calcualtions, but the static version can be used in its own right to calculate ground states, or other states of interest via suitable constraints.

The OneApi branch is created to make use of the Intel OneAPI tools which will allow GPU offloading including for the fftw and blas/lapack library routies.

