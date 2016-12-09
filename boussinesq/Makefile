F90 =  mpif90
#F90FLAGS =  -llapack -lblas -fconvert=big-endian -O3 -Wall 
F90FLAGS =  -llapack -lblas -fconvert=big-endian -O3 -Wall\
            -fbounds-check -finit-real=nan  -finit-integer=nan
CDFFLAGS = $(shell nf-config --fflags --flibs)
#CDFFLAGS = -I/usr/local/include -L/usr/local/lib -lnetcdf 
MPIFLAGS = 

F2PYFLAGS = -llapack -lblas
#F2PY_MPIFLAGS = -I/usr/include/mpi -lmpi -lmpi_f77 
F2PY_MPIFLAGS = $(shell mpif90 --showme:link)

all: check_density.f90 density.o timing_module.o
	 $(F90) check_density.f90 density.o timing_module.o $(F90FLAGS) $(CDFFLAGS) \
            -o check_density.x
clean :
	 rm -f *.o *.mod

density.o: density.f90
	$(F90) $(F90FLAGS) -c density.f90
timing_module.o: timing_module.f90
	$(F90) $(F90FLAGS) -c timing_module.f90



