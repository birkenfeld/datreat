
SHELL = /bin/sh
SRCF=$(wildcard *.f) 
SRCF90=integration.f90 rpa_laplace.f90 polynom.f90 lmfit.f90 $(wildcard *.f90) 


OBJF=$(SRCF:.f=.o)
OBJF90=$(SRCF90:.f90=.o)

OBJ=$(OBJF) $(OBJF90)
.SUFFIXES:
.SUFFIXES: .f90 .f .o 
# FFLAGS=-march=core2 -O3  -save -fPIC -mfpmath=sse -I../
FFLAGS= -O3  -save -fPIC -I../
THEOS := $(shell grep -h "^ *[Ff][Uu][Nn][Cc][Tt][Ii][Oo][nN] .*( *x *, *pa *, *.*, *.*, *.*, *.*)" *.f *.f90 | sed 's>[Ff]......[Nn] *\(.*\)( *x.*)>\1>'|sort )
DOUBLE := $(shell grep -h "^ *[Ff][Uu][Nn][Cc][Tt][Ii][Oo][nN] .*( *x *, *pa *, *.*, *.*, *.*, *.*)" *.f *.f90  | sed 's>[Ff]......[Nn] *\(.*\)( *x.*)>\1>'|sort |uniq -d)

FC=gfortran
#FFLAGS= -O3 -ffree-line-length-none -fno-automatic -fPIC -fno-backslash -Wunused-function -Wunused-variable -I../ -fopenmp
FFLAGS= -O3 -fcheck=bounds -ffree-line-length-none -fno-automatic -fPIC -fno-backslash -Wunused-function -Wunused-variable -fopenmp -I../ 
# FFLAGS= -O3 -ffree-line-length-none -fno-automatic -fPIC -fno-backslash -Wunused-function -Wunused-variable -I../ 


#FFLAGS=-march=core2 -O2 -ffree-line-length-none -fno-automatic -fPIC -fno-backslash -msse4a -mfpmath=sse -Wunused-function -Wunused-variable -I../

#FC=ifort

libtheos.a: init_theories.o sum_theories_up.i90 $(OBJ)
	ar cr $@ $(OBJ) 

init_theories.o:init_theories.f90 sum_theories_up.i90
	@echo " ---->>>>>>   recompile  init_theories"
	$(FC) -c $(FFLAGS) $< 

.f90.o:
	$(FC) -c $(FFLAGS) $? 

.f.o:
	$(FC) -c $(FFLAGS) $? 



ini_theories.i90 sum_theories_up.i90: $(SRCF) $(SRCF90)
	@echo "double names in THEOS: $(DOUBLE)"
	@echo "found  names in THEOS: $(THEOS)"
ifeq "$(strip $(DOUBLE))" ""                      # no double names
	@echo " .... found and included"
	./make_ini_include >  ini_theories.i90
	./make_sum_up_include > sum_theories_up.i90
	
else                                # double names
$(error  "!!!Fehler- in $(DOUBLE) : function name used twice change one of it  !!!!!!!!!!  (double function names are not allowed)  )
endif

clean:
	rm -f *.o *.i90 libtheos* *.mod

