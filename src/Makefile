.SUFFIXES:
.SUFFIXES: .f90 .f .o .F90 .c
SHELL = /bin/sh
LOCAL:=$(shell cd ..;/bin/pwd)

ADDLIBGKS=-lGKS

SRCSTHEO=  $(wildcard theos/*.f) $(wildcard theos/*.f90) $(wildcard theos/*.i90)

FC=gfortran
#FFLAGS= -O3 -ffree-line-length-none -fno-automatic -fPIC -fno-backslash -Wunused-function -Wunused-variable -fopenmp 
FFLAGS= -O3 -ffree-line-length-none -fno-automatic -fPIC -fno-backslash -fopenmp -fcheck=bounds \
            -fallow-argument-mismatch -fallow-invalid-boz 
#FFLAGS= -O3 -ffree-line-length-none -fno-automatic -fPIC -fno-backslash -fcheck=bounds 

CC=gcc 

LDFLAGS = -DPATH_TO_LOCAL_DATREAT=\"$(LOCAL)\"  -DPATH_TO_MAKROS=\"$(LOCAL)/makros/\"  -DOPEN_MANUAL=\"open\ $(LOCAL)/doc/NewManualBook/ManualDtr.pdf\"   \
-DDEFAULT_EDIT=\"vim\" -Wl,-rpath,/usr/local/gr/lib -lGR 
#-static-libgfortran

OBJS = datreat_nc.o dataconvolution4.o dtpur.o util_get_y_value.o \
 new_com.o os_utils.o utilux.o  \
datreat_nc_subroutines.o PhysicalConstantsPlus.o commons_mod.o unift.o \
 linenoise.o linenoise_f90.o datreat_grframework.o 

LIBS = theos/libtheos.a \
lib/dtpurcollection/libdtpurcollection.a lib/fftpack5/libfftpack.a \
lib/slatec/libslatec.a lib/sminpack/libsminpack.a

MODS =  os_utils.F90 new_com.F90 PhysicalConstantsPlus.f90 commons_mod.f90 unift.f90 datreat_grframework.f90 theos/lmfit.f90

CCS = linenoise.c linenoise_f90.c

.F90.o:
	$(FC) -c $(FFLAGS) $(LDFLAGS) $(MODS) $<

.f90.o:
	$(FC) -c $(FFLAGS) $(LDFLAGS) $(MODS) $< 

.f.o:
	$(FC) -c $(FFLAGS) $(LDFLAGS) $(MODS) $< 

.c.o:
	$(CC) -c $(CCS) $< 


datreat: $(OBJS) $(LIBS) $(OBJ_GR)
	$(FC) -o datreat $(FFLAGS) $(OBJS) $(LIBS) $(LDFLAGS) $(OBJ_GR)  -L/usr/local/gr/lib $(ADDLIBGR) $(ADDLIBGKS)  
	@/bin/mv datreat ../bin
	@/bin/ln -sf $(LOCAL)/bin/datreat ~/bin
	@echo "___________________________________________________________________"
	@echo " making datreat was SUCCESSFUL                                     "
#	@head CHANGES.log
	@echo "___________________________________________________________________"

# CPP:	
#	@cpp commons-cpp.h commons.h
#       @cpp commons-cpp.h -o commons.h     fuer MAC OS-X weg damit
intel:
	make datreat FC=ifort FFLAGS="-g -O2" \
	LDFLAGS="-Wl,-rpath,/usr/local/intel/fc/stable/lib"

lib/dtpurcollection/libdtpurcollection.a:
	make -C lib/dtpurcollection FC=$(FC) FFLAGS="$(FFLAGS)"
lib/fftpack5/libfftpack.a:
	make -C lib/fftpack5 FC=$(FC) FFLAGS="$(FFLAGS)"
lib/slatec/libslatec.a:
	make -C lib/slatec FC=$(FC) FFLAGS="$(FFLAGS)"
lib/sminpack/libsminpack.a:
	make -C lib/sminpack FC=$(FC) FFLAGS="$(FFLAGS)"

theos/libtheos.a:$(SRCSTHEO)
	@echo "" ------------------------------------------------
	@echo "       Starting with Theories"
	$(MAKE) -C theos
	#cd theos/ && $(MAKE)
	@echo "=======>>>>>>>>>>     "
	@echo "       back from Theories "
	@echo "" ------------------------------------------------

#	make -C theos FC=$(FC) FFLAGS="$(FFLAGS)"


	

clean:
	make -C lib/dtpurcollection clean
	make -C lib/fftpack5 clean
	make -C lib/slatec clean
	make -C lib/sminpack clean
	make -C theos clean
	rm -f *.o
	rm -f *.mod
	rm -f unused_theos/*.o
	


cleanall:
	make -C lib/dtpurcollection cleanall
	make -C lib/fftpack5 cleanall
	make -C lib/slatec cleanall
	make -C lib/sminpack cleanall
	make -C theos clean
	rm -f *.o
	rm -f *.mod
	rm -f unused_theos/*.o


 
?:help 
