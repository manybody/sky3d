Sky3d is a nuclear time-dependent Hartree-Fock code.
It represents wave functions in a symmetry-unrestricted three-dimensional cartesian box, and operates in static and time-dependent modes.  The static mode is mainly for the generation of initial conditions for time-depdendent calcualtions, but the static version can be used in its own right to calculate ground states, or other states of interest via suitable constraints.

This is the v1.2 of the code which include an external boost of multipole type where the user can provide custom input that decides the nature of the multipole (monopole, quadrupole, octupole, and so on) boost. The principal aim is to calculate the multipole strength function, which is the Fourier transform of the time-dependent expectation value of the multipole operator, which has the same form as the external boost. The proper unit conversion is done so that we can extract the exact unit of the thus calculated strength function, which is comparable to the available literature in the field. The boundary conditions are chosen such that a Woods-Saxon-like function cuts off the external field, driving it to zero at the boundary. 

# To compile the following modules are necessary 

## * GCC/10.2.0 OpenMPI/4.0.5 FFTW/3.3.8 LAPACK/3.9.1 OpenBLAS/0.3.12 *