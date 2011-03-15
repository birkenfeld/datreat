
.SUFFIXES:
.SUFFIXES: .f .o

SS= $(wildcard *.f)
OBJ = $(SS:.f=.o)

libfftpack.a: $(OBJ)
	ar cr $@ $(OBJ)

.f.o:
	$(FC) -c $(FFLAGS) $< 

clean:
	rm -f *.o
cleanall:
	rm -f *.o libfftpack*
