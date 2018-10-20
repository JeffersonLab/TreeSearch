#------------------------------------------------------------------------------
# Core library

SRC  = Tracker.cxx Plane.cxx Hit.cxx Hitpattern.cxx \
	Projection.cxx Pattern.cxx PatternTree.cxx PatternGenerator.cxx \
	TreeWalk.cxx Node.cxx Road.cxx

EXTRAHDR = Helper.h Types.h EProjType.h

CORE = TreeSearch
CORELIB  = lib$(CORE).so
COREDICT = $(CORE)Dict

LINKDEF = $(CORE)_LinkDef.h

#------------------------------------------------------------------------------
# MWDC library

MWDCSRC  = MWDC.cxx WirePlane.cxx WireHit.cxx HitpatternLR.cxx \
	ProjectionLR.cxx TimeToDistConv.cxx BigBite.cxx

MWDC  = TreeSearch-MWDC
MWDCLIB  = lib$(MWDC).so
MWDCDICT = $(MWDC)Dict

MWDCLINKDEF = $(MWDC)_LinkDef.h

#------------------------------------------------------------------------------
# GEM library

GEMSRC  = GEMTracker.cxx GEMPlane.cxx GEMHit.cxx

GEM  = TreeSearch-GEM
GEMLIB  = lib$(GEM).so
GEMDICT = $(GEM)Dict

GEMLINKDEF = $(GEM)_LinkDef.h

#------------------------------------------------------------------------------
# Compile debug version (for gdb)
#export DEBUG = 1
# Compile extra code for printing verbose messages (enabled with fDebug)
export VERBOSE = 1
# Compile extra diagnostic code (extra computations and global variables)
export TESTCODE = 1
# Compile support code for MC input data
export MCDATA = 1

#export I387MATH = 1
export EXTRAWARN = 1

# Architecture to compile for
MACHINE := $(shell uname -s)
ARCH    := linux
ifeq ($(MACHINE),Darwin)
  ARCH := macosx
endif

#------------------------------------------------------------------------------
# Directory locations. All we need to know is INCDIRS.
# INCDIRS lists the location(s) of the C++ Analyzer header (.h) files

ifndef ANALYZER
  $(error $$ANALYZER environment variable not defined)
endif

INCDIRS  = $(wildcard $(addprefix $(ANALYZER)/, include src hana_decode HallA Podd))

ifdef EVIO_INCDIR
  INCDIRS += ${EVIO_INCDIR}
else ifdef EVIO
  INCDIRS += ${EVIO}/include
endif

#------------------------------------------------------------------------------
# Do not change anything  below here unless you know what you are doing

ifeq ($(strip $(INCDIRS)),)
  $(error No Analyzer header files found. Check $$ANALYZER)
endif

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTBIN      := $(shell root-config --bindir)
CXX          := $(shell root-config --cxx)
LD           := $(shell root-config --ld)

PKGINCLUDES  = $(addprefix -I, $(INCDIRS) ) -I$(shell pwd)
INCLUDES     = -I$(shell root-config --incdir) $(PKGINCLUDES)

LIBS          = 
GLIBS         = 

ifeq ($(ARCH),linux)
# Linux with gcc (RedHat)
ifdef DEBUG
  CXXFLAGS    = -g -O0
  LDFLAGS     = -g -O0
  DEFINES     =
else
  CXXFLAGS    = -O2 -g #-march=pentium4
  LDFLAGS     = -O -g
  DEFINES     = -DNDEBUG
endif
DEFINES      += -DLINUXVERS -DHAS_SSTREAM
CXXFLAGS     += -Wall -Woverloaded-virtual -fPIC
DICTCXXFLG   :=
ifdef EXTRAWARN
#FIXME: should be configure'd:
CXXVER       := $(shell g++ --version | head -1 | sed 's/.* \([0-9]\)\..*/\1/')
ifeq ($(CXXVER),4)
CXXFLAGS     += -Wextra -Wno-missing-field-initializers -Wno-unused-parameter
DICTCXXFLG   := -Wno-strict-aliasing 
endif
endif
SOFLAGS       = -shared
ifdef I387MATH
CXXFLAGS     += -mfpmath=387
else
CXXFLAGS     += -march=core2 -mfpmath=sse
endif
endif

ifeq ($(ARCH),macosx)
# Mac OS X with gcc >= 3.x or clang++ >= 5
ifdef DEBUG
  CXXFLG     := -g -O0
  LDFLAGS    := -g -O0
  DEFINES    :=
else
  CXXFLG     := -O
  LDFLAGS    := -O
  DEFINES    := -DNDEBUG
endif
DEFINES      += -DMACVERS
CXXFLG       += -Wall -fPIC
CXXEXTFLG     = -Woverloaded-virtual
LD           := $(CXX)
LDCONFIG     :=
SOFLAGS      := -shared -Wl,-undefined,dynamic_lookup
SONAME       := -Wl,-install_name,
ifeq ($(CXX),clang++)
CXXEXTFLG    += -Wextra -Wno-missing-field-initializers -Wno-unused-parameter
else
#FIXME: should be configure'd:
CXXVER       := $(shell g++ --version | head -1 | sed 's/.* \([0-9]\)\..*/\1/')
ifeq ($(CXXVER),4)
CXXEXTFLG    += -Wextra -Wno-missing-field-initializers -Wno-unused-parameter
DICTCXXFLG   := -Wno-strict-aliasing
endif
endif
endif

ifdef VERBOSE
DEFINES      += -DVERBOSE
endif
ifdef TESTCODE
DEFINES      += -DTESTCODE
endif
ifdef MCDATA
DEFINES      += -DMCDATA
endif

CXXFLAGS     += $(DEFINES) $(ROOTCFLAGS) $(ROOTCFLAGS) $(PKGINCLUDES)
LIBS         += $(ROOTLIBS) $(SYSLIBS)
GLIBS        += $(ROOTGLIBS) $(SYSLIBS)

ifndef PKG
PKG           = lib$(CORE)
LOGMSG        = "$(PKG) source files"
else
LOGMSG        = "$(PKG) Software Development Kit"
endif
DISTFILE      = $(PKG).tar

#------------------------------------------------------------------------------
OBJ           = $(SRC:.cxx=.o) $(COREDICT).o
HDR           = $(SRC:.cxx=.h) $(EXTRAHDR)
DEP           = $(SRC:.cxx=.d)

MOBJ          = $(MWDCSRC:.cxx=.o) $(MWDCDICT).o
MHDR          = $(MWDCSRC:.cxx=.h)
MDEP          = $(MWDCSRC:.cxx=.d)

GOBJ          = $(GEMSRC:.cxx=.o) $(GEMDICT).o
GHDR          = $(GEMSRC:.cxx=.h)
GDEP          = $(GEMSRC:.cxx=.d)

all:		$(CORELIB) $(MWDCLIB) $(GEMLIB)

mwdc:		$(MWDCLIB)

gem:		$(GEMLIB)

$(CORELIB):	$(OBJ)
		$(LD) $(LDFLAGS) $(SOFLAGS) -o $@ $^
		@echo "$@ done"

$(MWDCLIB):	$(MOBJ) $(CORELIB)
		$(LD) $(LDFLAGS) $(SOFLAGS) -o $@ $^
		@echo "$@ done"

$(GEMLIB):	$(GOBJ) $(CORELIB)
		$(LD) $(LDFLAGS) $(SOFLAGS) -o $@ $^
		@echo "$@ done"

ifeq ($(ARCH),linux)
$(COREDICT).o:	$(COREDICT).cxx
	$(CXX) $(CXXFLAGS) $(DICTCXXFLG) -o $@ -c $^
$(MWDCDICT).o:	$(MWDCDICT).cxx
	$(CXX) $(CXXFLAGS) $(DICTCXXFLG) -o $@ -c $^
$(GEMDICT).o:	$(GEMDICT).cxx
	$(CXX) $(CXXFLAGS) $(DICTCXXFLG) -o $@ -c $^
endif

$(COREDICT).cxx: $(HDR) $(LINKDEF)
	@echo "Generating dictionary $(COREDICT)..."
	$(ROOTBIN)/rootcint -f $@ -c $(INCLUDES) $(DEFINES) $^

$(MWDCDICT).cxx: $(MHDR) $(MWDCLINKDEF)
	@echo "Generating dictionary $(MWDCDICT)..."
	$(ROOTBIN)/rootcint -f $@ -c $(INCLUDES) $(DEFINES) $^

$(GEMDICT).cxx: $(GHDR) $(GEMLINKDEF)
	@echo "Generating dictionary $(GEMDICT)..."
	$(ROOTBIN)/rootcint -f $@ -c $(INCLUDES) $(DEFINES) $^

install:	all
		$(error Please define install yourself)
# for example:
#		cp $(USERLIB) $(LIBDIR)

clean:
		rm -f *.o *~ $(CORELIB) $(COREDICT).*
		rm -f $(MWDCLIB) $(MWDCDICT).* $(GEMLIB) $(GEMDICT).*

realclean:	clean
		rm -f *.d

srcdist:
		rm -f $(DISTFILE).gz
		rm -rf $(PKG)
		mkdir $(PKG)
		cp -p $(SRC) $(HDR) $(LINKDEF) db*.dat Makefile $(PKG)
		cp -p $(MWDCLINKDEF) $(GEMLINKDEF) $(SOLIDLINKDEF) $(PKG)
		cp -p $(MWDCSRC) $(MHDR) $(GEMSRC) $(GHDR) $(PKG)
		gtar czvf $(DISTFILE) --ignore-failed-read \
		 -V $(LOGMSG)" `date -I`" $(PKG)
		rm -rf $(PKG)

develdist:	srcdist
		mkdir $(PKG)
		ln -s ../.git $(PKG)
		cp -p .gitignore $(PKG)
		gunzip -f $(DISTFILE).gz
		gtar rhvf $(DISTFILE) --exclude=*~ $(PKG)
		xz -f $(DISTFILE)
		rm -rf $(PKG)

.PHONY: all clean realclean srcdist

.SUFFIXES:
.SUFFIXES: .c .cc .cpp .cxx .C .o .d

%.o:	%.cxx
ifeq ($(strip $(MAKEDEPEND)),)
	$(CXX) $(CXXFLAGS) -MMD -o $@ -c $<
	@mv -f $*.d $*.d.tmp
else
	$(CXX) $(CXXFLAGS) -o $@ -c $<
	$(MAKEDEPEND) $(ROOTINC) $(INCLUDES) $(DEFINES) -c $< > $*.d.tmp
endif
	@sed -e 's|.*:|$*.o:|' < $*.d.tmp > $*.d
	@sed -e 's/.*://' -e 's/\\$$//' < $*.d.tmp | fmt -1 | \
	  sed -e 's/^ *//' -e 's/$$/:/' >> $*.d
	@rm -f $*.d.tmp

###

-include $(DEP)
-include $(MDEP)
-include $(GDEP)

