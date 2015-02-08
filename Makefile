.SUFFIXES: .f90 
obj =  module.o main.o readHmnR2.o \
       ham_qlayer2qlayer.o inverse.o \
		 chernrealspace.o chern.o ham_bulk.o mat_mul.o \
       readinput.o  ham_twist.o detmatrix.o ham_haldane.o \
       disorder.o eigen.o ek_bulk.o
cuda_obj = # cuda_face.o

FC  = mpif90 # -warn all -check all -traceback -pg
#FC  = ifort
CPPFLAGS= #-D__CUDA__  -D__CULA__ # -D__MAGMA__
FLAGS = -fpp -O3 -nogen-interface $(CPPFLAGS) # -check all

 LIB = -L. -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread
#LIB = -L. -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -liomp5 -lpthread


Jupiter :  $(obj) $(cuda_obj)
	$(FC) $(obj) $(cuda_obj) -o Jupiter $(LIB)

.f90.o: 
	$(FC) -c $(FLAGS) $*.f90 

clean :
	rm -f *.o *.mod *~ Jupiter
