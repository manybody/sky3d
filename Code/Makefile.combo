# Compiler
seq omp seq_debug omp_debug debug  : export COMPILER = gfortran
mpi mpi_debug mpi-omp              : export COMPILER = mpif90
ifort ifort_seq ifort_omp          : export COMPILER = ifort
cuda cuda_debug                    : export COMPILER = nvfortran

# Compiler flags
seq mpi                   : export FLAGS = -cpp -O3 -msse4.2 -mfpmath=sse -ffast-math -finline-functions -funroll-loops
omp mpi-omp               : export FLAGS = -cpp -O3 -fopenmp -msse4.2 -mfpmath=sse -ffast-math -finline-functions -funroll-loops
seq_debug mpi_debug debug : export FLAGS = -cpp -g -fbacktrace
omp_debug                 : export FLAGS = -cpp -g -fbacktrace -fopenmp
ifort ifort_seq           : export FLAGS = -cpp -axsse4.2 -msse4.2 -O3 -ip -no-prec-div -align all -mcmodel=large -shared-intel
ifort_omp                 : export FLAGS = -cpp -axsse4.2 -msse4.2 -O3 -ip -no-prec-div -align all -mcmodel=large -shared-intel -qopenmp
cuda                      : export FLAGS = -cpp -D CUDA -O3 -msse4.2 -Mcuda -mcmodel=medium -Mlarge_arrays -r8
cuda_debug                : export FLAGS = -cpp -D CUDA -g -traceback -O3 -msse4.2 -Mcuda -mcmodel=medium -Mlarge_arrays -r8

# Libraries
seq debug seq_debug mpi mpi_debug mpi_omp omp omp_debug ifort ifort_seq ifort_omp : export LIBS = -lfftw3 -llapack
cuda cuda_debug                                                                   : export LIBS = -lcufft -lcufftw -llapack

# Building object file
seq omp seq_debug omp_debug debug ifort ifort_seq ifort_omp cuda cuda_debug : export PARALLEL = sequential.o
mpi mpi_debug mpi-omp                                                       : export PARALLEL = parallel.o

# CUDA object files
cuda cuda_debug : export CUDA_OBJ = cufft_m.o precision_m.o

# Executables
seq debug seq_debug : export EXEC = sky3d.seq
ifort ifort_seq     : export EXEC = sky3d.ifort.seq
mpi mpi_debug       : export EXEC = sky3d.mpi
omp omp_debug       : export EXEC = sky3d.omp
ifort_omp           : export EXEC = sky3d.ifort.omp
mpi-omp             : export EXEC = sky3d.mpi-omp
cuda cuda_debug     : export EXEC = sky3d.cuda

# Object files
OBJS = params.o grids.o levels.o fourier.o forces.o $(PARALLEL) \
	inout.o coulomb.o trivial.o densities.o fragments.o twobody.o \
	energies.o static.o meanfield.o dynamic.o pairs.o moment.o \
	main3d.o user.o external.o abso_bc.o $(CUDA_OBJ)

SOURCE_TRAILER = .f90
SHELL = /bin/sh

# Clean
clean:
	-@rm -f *.o *.mod 2>/dev/null || true

clean-exec:
	-@rm -f *.o *.mod sky3d.* 2>/dev/null || true

# Build
all:$(EXEC)
$(EXEC):$(OBJS)
	$(COMPILER) $(FLAGS) $(LIBS) -o $@ $(OBJS)

user.o : user.f90 params.o grids.o levels.o

external.o : external.f90 params.o $(PARALLEL) grids.o levels.o meanfield.o

coulomb.o : coulomb.f90 params.o grids.o densities.o

inout.o : inout.f90 params.o $(PARALLEL) grids.o moment.o densities.o \
	levels.o  coulomb.o forces.o meanfield.o forces.o

trivial.o : trivial.f90 params.o grids.o

forces.o : forces.f90 forces.data params.o

fourier.o : fourier.f90 $(CUDA_OBJ) params.o grids.o

params.o : params.f90

moment.o : moment.f90 params.o grids.o densities.o

sequential.o: sequential.f90 levels.o params.o grids.o densities.o

parallel.o: parallel.f90 levels.o params.o grids.o densities.o

levels.o : levels.f90 params.o grids.o forces.o fourier.o trivial.o

grids.o : grids.f90 params.o forces.o

densities.o : densities.f90 params.o grids.o levels.o trivial.o

fragments.o : fragments.f90 params.o grids.o forces.o levels.o twobody.o \
	$(PARALLEL)

twobody.o : twobody.f90 params.o grids.o densities.o moment.o forces.o

energies.o :   energies.f90 params.o forces.o densities.o levels.o grids.o \
	trivial.o pairs.o

meanfield.o : meanfield.f90 params.o densities.o forces.o grids.o coulomb.o \
	trivial.o levels.o

pairs.o : pairs.f90 forces.o grids.o levels.o densities.o

static.o : static.f90 params.o densities.o grids.o \
	levels.o moment.o energies.o inout.o pairs.o meanfield.o

dynamic.o : dynamic.f90 params.o densities.o grids.o trivial.o \
	levels.o moment.o energies.o inout.o meanfield.o $(PARALLEL) \
	twobody.o external.o abso_bc.o

main3d.o : main3d.f90 params.o fourier.o forces.o densities.o meanfield.o \
    levels.o grids.o fragments.o $(PARALLEL) dynamic.o static.o coulomb.o user.o

.SUFFIXES:  .f90 .f .o

%.o : %.mod

%.o : %.f90
	$(COMPILER) $(FLAGS) -c $<

%.o : %f
	$(COMPILER) $(FLAGS) -c $<
