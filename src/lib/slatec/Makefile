
.SUFFIXES:
.SUFFIXES: .f .o

#OBJ =  aaaaaa.o \
#	acosh.o \
#	aie.o \
#	ai.o \
#	spofa.o \
#	spofs.o \
#	spoir.o \

SS= $(wildcard *.f)
OBJ = $(SS:.f=.o)

libslatec.a: $(OBJ)
	ar cr $@ $(OBJ)

.f.o:
	$(FC) -c $(FFLAGS) $< 

clean:
	rm -f *.o
cleanall:
	rm -f *.o libslatec*
