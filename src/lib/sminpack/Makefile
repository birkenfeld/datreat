
.SUFFIXES:
.SUFFIXES: .f .o

#OBJ  =   enorm.o fdjac2.o lmdif1.o lmdif.o lmpar.o qrfac.o qrsolv.o spmpar.o
SS= $(wildcard *.f)
OBJ = $(SS:.f=.o)
libsminpack.a: $(OBJ)
	ar cr $@ $(OBJ)

.f.o:
	$(FC) -c $(FFLAGS) $< 

clean:
	rm -f *.o
cleanall:
	rm -f *.o libsminpack*
