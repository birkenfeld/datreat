
.SUFFIXES:
.SUFFIXES: .f .o

#OBJ  =   642_TOMS_spline.o s09hre.o sadhre.o sfshre.o sgetri.o srlhre.o \
# s113re.o schhre.o sgetf2.o sinhre.o strhre.o s07hre.o s132re.o scuhre.o \
# sgetrf.o slaswp.o strtri.o

SS= $(wildcard *.f)
OBJ = $(SS:.f=.o)

libdtpurcollection.a: $(OBJ)
	ar cr $@ $(OBJ)

.o.o:
	$(FC) -c $(FFLAGS) $< 

cleanall:
	rm -f *.o libdtpurcollection*
clean:
	rm -f *.o 
