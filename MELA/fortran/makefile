ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
	false
else

RM  = /bin/rm

MAINDIR = $(shell pwd)/
MELALIBDIR = ${MELA_LIB_PATH}
OBJECTSDIR=
LIB = libjhugenmela.so
LIBRULE = $(MELALIBDIR)/$(LIB)

INCLUDE = -I$(MAINDIR)
Comp = gfort

# Linking LHAPDF
# link pdfs via LHA library ('Yes' or 'No')
UseLHAPDF=No
# directory which contains libLHAPDF.a, libLHAPDF.la, libLHAPDF.so
LHAPDFLIBDIR=/./
# remember to export
#          LD_LIBRARY_PATH=/.../LHAPDF-x.y.z/lib/:${LD_LIBRARY_PATH}
#          LHAPDF_DATA_PATH=/.../LHAPDF-x.y.z/share/LHAPDF/:${LHAPDF_DATA_PATH}
ifeq ($(UseLHAPDF),Yes)
   LHAPDFflags = -L$(LHAPDFLIBDIR) -lLHAPDF -DuseLHAPDF=1
else
   LHAPDFflags = -DuseLHAPDF=0
endif

# Linking the Collier library
UseCOLLIER=Yes
# directory which contains libLHAPDF.a, libLHAPDF.la, libLHAPDF.so
COLLIERLIBDIR=$(MAINDIR)../data/${MELA_ARCH}/
COLLIERINCDIR=$(MAINDIR)../COLLIER/
# remember to export
#          LD_LIBRARY_PATH=/.../LHAPDF-x.y.z/lib/:${LD_LIBRARY_PATH}
#          LHAPDF_DATA_PATH=/.../LHAPDF-x.y.z/share/LHAPDF/:${LHAPDF_DATA_PATH}
ifeq ($(UseCOLLIER),Yes)
   COLLIERflags = -L$(COLLIERLIBDIR) -lcollier -DuseCollier=1 -I$(COLLIERINCDIR)
else
   COLLIERflags = -DuseCollier=0
endif


ifeq ($(Comp),ifort)
	fcomp = ifort -fpp -O2 -vec-report0 -Dcompiler=1 -lifcore $(LHAPDFflags) $(COLLIERflags) -fPIC
endif
ifeq ($(Comp),gfort)
	fcomp = gfortran -O0 -ffree-line-length-none -Dcompiler=2 $(LHAPDFflags) $(COLLIERflags) -fno-automatic -fno-f2c  -fPIC -g 
endif

ifeq ($(Comp),gfort)
	ccomp = gcc
endif
ifeq ($(Comp),ifort)
	ccomp = icc
endif

ifeq ($(Comp),gfort)
	clib =  -lm -lgfortran
endif
ifeq ($(Comp),ifort)
	clib =  -lm -lirc
endif


all:    $(LIBRULE)

.SILENT:


SOURCESMAINF90 = $(wildcard $(MAINDIR)*.F90)
SOURCESMAINFF = $(wildcard $(MAINDIR)*.f)
SOURCESCPS = $(wildcard $(MAINDIR)CPS/*.f)
SOURCESPDF = $(wildcard $(MAINDIR)pdfs/*.f)

OBJECTSMAINF90PRIM = $(SOURCESMAINF90:.F90=.o)
OBJECTSMAINFFPRIM = $(SOURCESMAINFF:.f=.o)
OBJECTSCPSPRIM = $(SOURCESCPS:.f=.o)
OBJECTSPDFPRIM = $(SOURCESPDF:.f=.o)

OBJECTSMAINF90 = $(subst $(MAINDIR),$(OBJECTSDIR),$(OBJECTSMAINF90PRIM))
OBJECTSMAINFF = $(subst $(MAINDIR),$(OBJECTSDIR),$(OBJECTSMAINFFPRIM))
OBJECTSCPS = $(subst $(MAINDIR)CPS/,$(OBJECTSDIR),$(OBJECTSCPSPRIM))
OBJECTSPDF = $(subst $(MAINDIR)pdfs/,$(OBJECTSDIR),$(OBJECTSPDFPRIM))

OBJECTS = $(OBJECTSCPS) $(OBJECTSPDF) $(OBJECTSMAINFF) $(OBJECTSMAINF90)


# General rules for object files
$(OBJECTSDIR)%.o: %.F90 | $(OBJECTSDIR)mod_Parameters.o $(OBJECTSDIR)mod_Misc.o
	@echo " Compiling $< with $(Comp)"; \
	$(fcomp) -c $<
$(OBJECTSDIR)%.o: %.f
	@echo " Compiling $< with $(Comp)"; \
	$(fcomp) -c $<
$(OBJECTSDIR)%.o: CPS/%.f
	@echo " Compiling $< with $(Comp)"; \
	$(fcomp) -c $<
$(OBJECTSDIR)%.o: pdfs/%.f
	@echo " Compiling $< with $(Comp)"; \
	$(fcomp) -c $<

# Special rules for modules with different dependencies
$(OBJECTSDIR)mod_Parameters.o: mod_Parameters.F90
	@echo " Compiling $< with $(Comp)"; \
	$(fcomp) -c $<
$(OBJECTSDIR)mod_Misc.o: mod_Misc.F90 $(OBJECTSDIR)mod_Parameters.o
	@echo " Compiling $< with $(Comp)"; \
	$(fcomp) -c $<
$(OBJECTSDIR)mod_TTBHiggs.o: mod_TTBHiggs.F90 $(OBJECTSDIR)mod_Parameters.o $(OBJECTSDIR)mod_Misc.o $(OBJECTSDIR)mod_TopDecay.o
	@echo " Compiling $< with $(Comp)"; \
	$(fcomp) -c $<
$(OBJECTSDIR)mod_THiggs.o: mod_THiggs.F90 $(OBJECTSDIR)mod_Parameters.o $(OBJECTSDIR)mod_Misc.o $(OBJECTSDIR)mod_TopDecay.o
	@echo " Compiling $< with $(Comp)"; \
	$(fcomp) -c $<
$(OBJECTSDIR)mod_JHUGenMELA.o: mod_JHUGenMELA.F90 $(OBJECTSDIR)mod_Parameters.o $(OBJECTSDIR)mod_Misc.o $(OBJECTSDIR)mod_Kinematics.o
	@echo " Compiling $< with $(Comp)"; \
	$(fcomp) -c $<
$(OBJECTSDIR)mod_JHUGen.o: mod_JHUGen.F90 $(OBJECTSDIR)mod_HashCollection.o $(OBJECTSDIR)mod_JHUGenMELA.o
	@echo " Compiling $< with $(Comp)"; \
	$(fcomp) -c $<

$(LIBRULE): $(OBJECTS)
	@echo "Linking $(LIB)"; \
        g++ -Wl,-rpath=$(MELALIBDIR),-soname,$(LIB) -L$(MELALIBDIR) -shared *.o -o $@

clean:
	rm -f *.so $(LIBRULE) $(OBJECTSDIR)*.o $(OBJECTSDIR)*.mod


endif
