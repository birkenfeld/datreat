
.SUFFIXES:
.SUFFIXES: .f90 .f .o .F90
SHELL = /bin/sh
LOCAL:=$(shell cd ..;/bin/pwd)
##  To include GR and GKS if installed
ifneq "$(wildcard /usr/local/lib/libgks.*)" ""
      ADDLIBGKS=-lgks
else ifneq "$(wildcard /usr/local/gr/lib/libGKS.*)" ""
      ADDLIBGKS=-lgks
else
      ADDLIBGKS=
endif
ifneq "$(wildcard /usr/local/lib/libgr.*)" ""
      ADDLIBGR=-lgr
	  OBJ_GR=datreat_plot2gr.o
else ifneq "$(wildcard /usr/local/gr/lib/libgr*.*)" ""
      ADDLIBGR=-lgrsoft_ifort
	  OBJ_GR=datreat_plot2gr.o	  
else
      ADDLIBGR=
	  OBJ_GR=datreat_plot2gr_dummy.o
endif
#
# 
#
#  path to XmGrace (to libgrace_np.a)  and takes the first found 
#  standard path is  /usr/local/grace/lib  or in ubuntu /usr/lib/)
#  ---- add your path at the end if dirs here
GRACEdirs=/usr/local/grace/lib /usr/lib
find_files = $(wildcard $(dir)/libgrace_np.a)
GRACE= $(dir $(firstword $(foreach dir,$(GRACEdirs),$(find_files))   ))
#-------------------------------------------
#due to a bug in grace we need to use the old grace_np library of version 5.1.20
#-L/usr/local/grace/lib -lgrace_np-5.1.20
#instead of the new one
#-lgrace_np
# so if the old one is there (iff filesystem) use this instead of the old
GRACENP=grace_np
ifeq ($(notdir $(wildcard $(GRACE)/libgrace_np-5.1.20.*)),libgrace_np-5.1.20.a)
	GRACENP=grace_np-5.1.20
endif

ifeq "$(GRACE)" ""   
   $(warning "                 READ this ")
   $(warning "++++++++++++++++++READ this+++++++++++++++++++")
   $(warning +$(shell  ls -l `whereis xmgrace | cut -d " " -f 2 `)) 
   $(warning "libgrace_np.a  not found in GRACE=/path/to/grace/ ,  Please add path to libgrace_np.a to GRACEdirs in  Makefile )
   $(warning " We need the path to the installation directory where  you find lib/libgrace_np , so follow the link given above ")
   $(warning "++++++++++++++++++READ this+++++++++++++++++++")
   $(warning "                    READ this             ")
   $(error   )
endif

SRCSTHEO=  $(wildcard theos/*.f) $(wildcard theos/*.f90) $(wildcard theos/*.i90)

#FC=gfortran
#FFLAGS=-march=core2 -O2 -ffree-line-length-none -fno-automatic -fPIC -fno-backslash -msse4a -mfpmath=sse -I$LOCAL -Wunused-function -Wunused-variable
FC=ifort
FFLAGS=-march=core2 -O2  -save -fPIC  -mfpmath=sse -I$LOCAL -heap_arrays 100 

#-fbounds-check
LDFLAGS = -DPATH_TO_LOCAL_DATREAT=\"$(LOCAL)\"


OBJS = datreat05_2.o dataconvolution4.o dtpur.o util_get_y_value.o \
 incom.o utilux.o datreat2xmgrace.o incom_subroutines.o \
datreat_subroutines.o PhysicalConstantsPlus.o commons_mod.o unift.o 

LIBS = theos/libtheos.a \
lib/dtpurcollection/libdtpurcollection.a lib/fftpack5/libfftpack.a \
lib/slatec/libslatec.a lib/sminpack/libsminpack.a

MODS = PhysicalConstantsPlus.f90 commons_mod.f90 unift.f90

.F90.o:
	$(FC) -c $(FFLAGS) $(LDFLAGS) $(MODS) $<

.f90.o:
	$(FC) -c $(FFLAGS) $(LDFLAGS) $(MODS) $< 

.f.o:
	$(FC) -c $(FFLAGS) $(LDFLAGS) $(MODS) $< 


datreat: $(OBJS) $(LIBS) $(OBJ_GR)
	$(FC) -o datreat $(FFLAGS) $(OBJS) $(LIBS) $(LDFLAGS) -L$(GRACE) $(OBJ_GR) -l$(GRACENP)  -L/usr/local/gr/lib $(ADDLIBGR) $(ADDLIBGKS) -Wl,-rpath,/usr/local/intel/lib/intel64 -lifcore -limf -lintlc -llapack_ifort -lblas_ifort  
	@/bin/mv datreat ../bin
	@/bin/ln -sf $(LOCAL)/bin/datreat ~/bin
	@echo "___________________________________________________________________"
	@echo " If ~/bin exists and is in your PATH start it from anywhere        "
	@echo "      clean local directory by >   make clean                      "
	@echo "___________________________________________________________________"
	@echo "Start with > datreat        "
	@echo "___________________________________________________________________"
	@echo "START 'example' to test everything_________________________________"
	@echo "___________________________________________________________________"
	@echo "___________________________________________________________________"
	@head CHANGES.log
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



distribution:	
	cd ../..; tar --exclude=*.o --exclude=*.a --exclude=*.so --exclude=*~ --exclude=bin/datreat --exclude=.git* -z -cf ndatreat_$(shell date +%F).tar.gz ndatreat ;
	cd ../..; ln -s ndatreat_$(shell date +%F).tar.gz ndatreat_latest.tar.gz
	touch ./datreat05_2.f90
	@if [ $$USERNAME = biehl  ]; then  mkdir -p /temp/biehl; cp ../../ndatreat_$(shell date +%F).tar.gz /temp/biehl/;chmod ugo+r /temp/biehl/ndatreat_$(shell date +%F).tar.gz;fi
	#chmod ugo+r /temp/biehl/ndatreat_$(shell date +%F).tar.gz

grace:
	@echo "to get gracerc files for startup do "
	@echo mkdir ~/.grace
	@echo cp -r ../grace/* ~/.grace

help: 
	@echo "make <Options>"
	@echo "Options:"
	@echo "            : builds datreat "
	@echo "clean       : delets objekt files executables were not touched"
	@echo "cleanall    : deletes evrything exept source"
	@echo "help        : shows this "
	@echo "distribution   : make a tar.gz to distribute to others with pure sourcefiles "
	@echo "------------ datreat by M. Monkenbusch O. Holderer R. Biehl"

?:help
