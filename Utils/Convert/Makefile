# for Intel
COMPILER = gfortran -frecord-marker=4
LINKER = gfortran -frecord-marker=4
#for debug
#COMPILERFLAGS  =    -g
#LINKERFLAGS    =  -g 
#for pentium 4
#COMPILERFLAGS  =    -O3 -fno-alias -xW -tpp7
#LINKERFLAGS    =  -O3  -fno-alias -xW -tpp7
#for pentium 3
COMPILERFLAGS  =    -O 
LINKERFLAGS    =  -O 

SOURCE_TRAILER = .f90
EXEC           = Tdhf2Silo
OBJS           = Tdhf2Silo.o
STDLIB = -lsilo
SHELL = /bin/sh

#---------------------------------------------------------	
all:$(EXEC) 
$(EXEC):$(OBJS) 	
	$(LINKER) $(LINKERFLAGS) -o $@  $(OBJS) $(STDLIB)
#
Tdhf2Silo.o : Tdhf2Silo.f90

.SUFFIXES:  .f90 .f .o 

.f90.o:	
	$(COMPILER) $(COMPILERFLAGS)  -c $<

.f.o:
	$(COMPILER) $(COMPILERFLAGS)  -c $<
