# All targets with # symbol are self-documenting, i.e. make help or simply make will
# show the targets among available options
#
# User targets are at the bottom
#
ifndef ROOTSYS
all:
	@echo "ROOTSYS is not set. Please set ROOT environment properly"; echo
else
ifndef ROOFITSYS
all:
	@echo "ROOFITSYS is not set. Please set ROOT environment properly"; echo
else


RM  = /bin/rm

MELADIR = $(shell pwd)
MELASRCDIR = $(MELADIR)/src
MELAOBJDIR = $(MELADIR)/obj
# Modify MELALIBDIR for the gcc version as needed
MELALIBDIR = ${MELA_LIB_PATH}
MADGRAPHLIBDIR = $(MELADIR)/data/MADGRAPH_PROCESSES
# _melapkgpath_ should refer to the root compilation path just like MELADIR with an extra '/'.
# If environment variables need to be inserted without expansion for portability,
# you can use '.oODOLLAROo..oOOPEN_BRACKETOo.[YOUR_ENV_VARIABLE].oOCLOSE_BRACKETOo.',
# e.g. '.oODOLLAROo..oOOPEN_BRACKETOo.CMSSW_BASE.oOCLOSE_BRACKETOo.' in place of '${CMSSW_BASE}' as used in the CMS-specific compilation.
MELAPKGPATH = $(MELADIR)/
LIB = libJHUGenMELAMELA.so
LIBRULE = $(MELALIBDIR)/$(LIB)

pyLIB = Mela$(shell python3-config --extension-suffix)
pyLIBRULE = $(MELALIBDIR)/$(pyLIB)
pyLINKFILE = $(MELADIR)/python/mela_binding.cpp

ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs) -Lrootlib
NLIBS = $(ROOTLIBS)
NLIBS += -L$(ROOFITSYS)/lib -lMinuit -lRooFitCore -lRooFit -lgfortran
LIBS = $(filter-out -lNew, $(NLIBS))

ROOFITINC = -I$(ROOFITSYS)/include
ROOTINC = $(ROOFITINC)

CC = g++
CPPINC = -I$(MELADIR)/interface $(ROOTINC)
CPPLOAD = -L$(MELALIBDIR) -lmcfm_709 -ljhugenmela -lcollier -L$(MADGRAPHLIBDIR) -lggF_SMEFTsim $(LIBS)
CPPOPTS =  -fPIC -O2 -ftree-vectorize -fipa-pta -felide-constructors -fvisibility-inlines-hidden -fno-math-errno \
	--param vect-max-version-for-alias-checks=50 -fmessage-length=0 -fdiagnostics-show-option \
	-Werror=unused-value \
	 $(shell root-config --cflags) -D_melapkgpath_=$(MELAPKGPATH)
CPPFLAGS = $(CPPOPTS) $(CPPINC) $(CPPLOAD)

LINKER = g++
LINKERFLAGS = -Wl,-rpath=$(MELALIBDIR),-soname,$(LIB) $(CPPLOAD)

SOURCESCC = $(wildcard $(MELASRCDIR)/*.cc)
SOURCESCXX = $(wildcard $(MELASRCDIR)/*.cxx)
OBJECTSPRIM = $(SOURCESCC:.cc=.o) $(SOURCESCXX:.cxx=.o)
OBJECTS = $(subst $(MELASRCDIR),$(MELAOBJDIR),$(OBJECTSPRIM))
DEPS = $(OBJECTS:.o=.d)

.PHONY: all help compile clean

.SILENT: alldirs $(OBJECTS) $(DEPS) $(MELAOBJDIR)/LinkDef_out.o $(LIBRULE) python/__init__.py $(pyLIBRULE) clean

all: $(LIBRULE) python/__init__.py $(pyLIBRULE)

nopython: $(LIBRULE) mad_test

$(MELAOBJDIR)/%.d:	$(MELASRCDIR)/%.c* | alldirs
	echo "Checking dependencies for $<"; \
	$(CC) -MM -MT $@ $(CPPFLAGS) $< > $@; \
                     [ -s $@ ] || rm -f $@

$(MELAOBJDIR)/%.o: 	$(MELASRCDIR)/%.c* $(MELAOBJDIR)/%.d | alldirs
	echo "Compiling $<"; \
	$(CC) $(CPPFLAGS) $< -c -o $@

$(MELAOBJDIR)/LinkDef_out.o: | alldirs
	echo "Making CINT dictionaries"; \
	rootcint -f $(MELASRCDIR)/LinkDef_out.cc -c -p $(CPPINC) $(MELASRCDIR)/LinkDef.h; \
	cat $(MELASRCDIR)/LinkDef.h $(MELASRCDIR)/LinkDef_out.cc > $(MELASRCDIR)/LinkDef_out.cxx; \
	rm $(MELASRCDIR)/LinkDef_out.cc; mv $(MELASRCDIR)/LinkDef_out_rdict.pcm $(MELALIBDIR)/LinkDef_out_rdict.pcm; \
	$(CC) $(CPPOPTS) $(CPPINC) -c $(MELASRCDIR)/LinkDef_out.cxx -o $(MELAOBJDIR)/LinkDef_out.o; \
	rm -f $(MELASRCDIR)/LinkDef_out.cxx

$(LIBRULE):	$(OBJECTS) $(MELAOBJDIR)/LinkDef_out.o | alldirs
	echo "Linking $(LIB)"; \
	$(LINKER) $(LINKERFLAGS) -shared $(OBJECTS) $(MELAOBJDIR)/LinkDef_out.o -o $@

$(pyLIBRULE): $(pyLINKFILE) $(OBJECTS) $(MELAOBJDIR)/LinkDef_out.o | $(LIBRULE)
	@echo "Making python"; \
	$(LINKER) $(pyLINKFILE) $(CPPFLAGS) $(LINKERFLAGS) \
	-O3 -Wall -shared $(OBJECTS) $(MELAOBJDIR)/LinkDef_out.o -std=c++17 -fPIC \
	$(shell python3 -m pybind11 --includes) -o $@

mad_test: test/testME_madgraph.cpp | alldirs
	@$(LINKER) $(CPPFLAGS) -L$(MELALIBDIR) -lJHUGenMELAMELA test/testME_madgraph.cpp -g -o testME_madgraph

alldirs:
	mkdir -p $(MELAOBJDIR)

python/__init__.py:
	touch python/__init__.py

clean:
	rm -rf $(MELAOBJDIR) $(MELASRCDIR)/LinkDef_out* $(MELALIBDIR)/$(LIB) $(MELALIBDIR)/LinkDef_out_rdict.pcm


include $(DEPS)


endif
endif
