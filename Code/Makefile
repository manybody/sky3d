#
# ----- targets--------------------------------------------------------------
#
.PHONY : seq mpi omp debug seq_debug omp_debug mpi_debug mpi-omp tabc tabc_seq tabc_omp ifort ifort_seq ifort_omp makeit clean clean-exec
#
# ----- define the compiler being used---------------------------------------
#
seq omp seq_debug omp_debug debug            : export COMPILER_SKY = gfortran
mpi mpi_debug tabc tabc_seq tabc_omp mpi-omp : export COMPILER_SKY = mpif90
ifort ifort_seq ifort_omp                    : export COMPILER_SKY = ifort
#
# ----- define libraries-----------------------------------------------------
#
seq omp debug seq_debug omp_debug tabc tabc_seq tabc_omp ifort ifort_seq ifort_omp : export LIBS_SKY = -lfftw3 -llapack -lblas
mpi mpi_debug mpi-omp                                                              : export LIBS_SKY = -lfftw3 -lscalapack -llapack -lblas
#
# ----- define compiler flags and linker flags-------------------------------
#
seq mpi tabc tabc_seq     : export COMPILERFLAGS_SKY = -O3 -msse4.2 -mfpmath=sse -ffast-math -finline-functions -funroll-loops
omp tabc_omp mpi-omp      : export COMPILERFLAGS_SKY = -O3 -fopenmp -msse4.2 -mfpmath=sse -ffast-math -finline-functions -funroll-loops
seq_debug mpi_debug debug : export COMPILERFLAGS_SKY = -g -fbacktrace
omp_debug                 : export COMPILERFLAGS_SKY = -g -fbacktrace -fopenmp
ifort ifort_seq           : export COMPILERFLAGS_SKY = -axsse4.2 -msse4.2 -O3 -ip -no-prec-div -align all -mcmodel=large -shared-intel
ifort_omp                 : export COMPILERFLAGS_SKY = -axsse4.2 -msse4.2 -O3 -ip -no-prec-div -align all -mcmodel=large -shared-intel -openmp
#
# ----- which communication object file to use for building------------------
#
seq omp seq_debug omp_debug debug ifort ifort_seq ifort_omp : export PARALLEL_SKY = sequential.o
mpi mpi_debug mpi-omp                                       : export PARALLEL_SKY = parallel.o
tabc tabc_seq tabc_omp                                      : export PARALLEL_SKY = tabc.o
#
# ----- which linear algebra object file to use for building-----------------
#
seq omp seq_debug omp_debug debug ifort ifort_seq ifort_omp tabc tabc_seq tabc_omp : export LINALG_SKY = linalg_lapack.o
mpi mpi_debug mpi-omp                                                              : export LINALG_SKY = linalg_scalapack.o
#
# ----- executables----------------------------------------------------------
#
seq debug seq_debug : export EXEC_SKY = sky3d.seq
ifort ifort_seq     : export EXEC_SKY = sky3d.ifort.seq
mpi mpi_debug       : export EXEC_SKY = sky3d.mpi
omp omp_debug       : export EXEC_SKY = sky3d.omp
ifort_omp           : export EXEC_SKY = sky3d.ifort.omp
mpi-omp             : export EXEC_SKY = sky3d.mpi-omp
tabc tabc_seq       : export EXEC_SKY = sky3d.tabc
tabc_omp            : export EXEC_SKY = sky3d.tabc_omp   
#
# ----- restart make or clean------------------------------------------------
#
seq omp mpi debug seq_debug mpi_debug omp_debug mpi-omp tabc tabc_seq tabc_omp ifort ifort_seq ifort_omp:
	@$(MAKE) -f $(firstword $(MAKEFILE_LIST)) makeit

clean:
	rm -f *.o *.mod

clean-exec:
	rm -fr *.o *.mod sky3d.*
#
#
# ***** after restart********************************************************
#
# ----- object files---------------------------------------------------------
#
OBJS   = params.o grids.o levels.o fourier.o forces.o $(PARALLEL_SKY) \
	inout.o coulomb.o trivial.o densities.o fragments.o twobody.o \
	energies.o static.o meanfield.o dynamic.o pairs.o moment.o \
	main3d.o user.o external.o abso_bc.o constraint.o $(LINALG_SKY) \
        radius.o temperature.o
#
LINKERFLAGS_SKY    =  $(COMPILERFLAGS_SKY)
#
SOURCE_TRAILER = .f90
#
SHELL = /bin/sh
#
#-----------------------------------------------------------------------------
makeit : $(EXEC_SKY)

$(EXEC_SKY):$(OBJS)
	$(COMPILER_SKY) $(LINKERFLAGS_SKY) -o $@  $(OBJS) $(LIBS_SKY) 
#

temperature.o : temperature.f90 levels.o

user.o : user.f90 params.o grids.o levels.o Makefile

linalg_lapack.o : linalg_lapack.f90 params.o grids.o levels.o $(PARALLEL_SKY) Makefile

linalg_scalapack.o : linalg_scalapack.f90 params.o grids.o levels.o $(PARALLEL_SKY) Makefile

radius.o : radius.f90 params.o grids.o levels.o densities.o trivial.o Makefile

external.o : external.f90 params.o $(PARALLEL_SKY) grids.o levels.o meanfield.o Makefile

coulomb.o : coulomb.f90 params.o grids.o densities.o Makefile

inout.o : inout.f90 params.o $(PARALLEL_SKY) grids.o moment.o densities.o \
	levels.o  coulomb.o forces.o meanfield.o forces.o Makefile

trivial.o : trivial.f90 params.o grids.o Makefile

forces.o : forces.f90 forces.data params.o  Makefile

fourier.o : fourier.f90 params.o grids.o Makefile

params.o : params.f90 Makefile

moment.o : moment.f90 params.o grids.o densities.o Makefile

sequential.o: sequential.f90 levels.o params.o grids.o densities.o Makefile

parallel.o: parallel.f90 levels.o params.o grids.o densities.o Makefile

levels.o : levels.f90 params.o grids.o forces.o fourier.o trivial.o Makefile

grids.o : grids.f90 params.o forces.o Makefile

densities.o : densities.f90 params.o grids.o levels.o trivial.o  Makefile

fragments.o : fragments.f90 params.o grids.o forces.o levels.o twobody.o \
	$(PARALLEL_SKY) Makefile

abso_bc.o : abso_bc.f90 params.o grids.o forces.o levels.o meanfield.o \
	fragments.o $(PARALLEL_SKY)  Makefile

twobody.o : twobody.f90 params.o grids.o densities.o moment.o forces.o Makefile

energies.o :   energies.f90 params.o forces.o densities.o levels.o grids.o \
	trivial.o pairs.o Makefile

meanfield.o : meanfield.f90 params.o densities.o forces.o grids.o coulomb.o \
	trivial.o levels.o Makefile

pairs.o : pairs.f90 forces.o grids.o levels.o densities.o $(PARALLEL_SKY) Makefile

static.o : static.f90 params.o densities.o grids.o radius.o \
	levels.o moment.o energies.o inout.o pairs.o \
        constraint.o meanfield.o $(LINALG_SKY) temperature.o Makefile

constraint.o : constraint.f90 params.o densities.o grids.o \
	levels.o Makefile

dynamic.o : dynamic.f90 params.o densities.o grids.o trivial.o \
	levels.o moment.o energies.o inout.o meanfield.o $(PARALLEL_SKY) \
	twobody.o external.o abso_bc.o Makefile

main3d.o : main3d.f90 params.o fourier.o forces.o densities.o meanfield.o \
    levels.o grids.o fragments.o $(PARALLEL_SKY) dynamic.o static.o coulomb.o user.o \
    Makefile

.SUFFIXES:  .f90 .f .o 

%.o : %.mod

%.o : %.f90
	$(COMPILER_SKY) $(COMPILERFLAGS_SKY) -c $<

%.o : %f
	$(COMPILER_SKY) $(COMPILERFLAGS_SKY) -c $<
