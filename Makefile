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
	ProjectionLR.cxx TimeToDistConv.cxx

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
export DEBUG = 1
#export PROFILE = 1
# Compile extra code for printing verbose messages (enabled with fDebug)
export VERBOSE = 1
# Compile extra diagnostic code (extra computations and global variables)
export TESTCODE = 1

export I387MATH = 1
export EXTRAWARN = 1

# Architecture to compile for
ARCH          = linux
#ARCH          = solarisCC5

#------------------------------------------------------------------------------
# Directory locations. All we need to know is INCDIRS.
# INCDIRS lists the location(s) of the C++ Analyzer header (.h) files

ifndef ANALYZER
  $(error $$ANALYZER environment variable not defined)
endif

INCDIRS  = $(wildcard $(addprefix $(ANALYZER)/, include src hana_decode hana_scaler))

#------------------------------------------------------------------------------
# Do not change anything  below here unless you know what you are doing

ifeq ($(strip $(INCDIRS)),)
  $(error No Analyzer header files found. Check $$ANALYZER)
endif

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)
ROOTBIN      := $(shell root-config --bindir)

INCLUDES      = $(ROOTCFLAGS) $(addprefix -I, $(INCDIRS) ) -I$(shell pwd)

LIBS          = 
GLIBS         = 

ifeq ($(ARCH),linux)
# Linux with gcc (RedHat)
CXX           = g++
ifdef DEBUG
  CXXFLAGS    = -g -O0
  LDFLAGS     = -g -O0
  DEFINES     =
else
  CXXFLAGS    = -O2 -march=pentium4
  LDFLAGS     = -O
  DEFINES     = -DNDEBUG
endif
DEFINES      += -DLINUXVERS -DHAS_SSTREAM
CXXFLAGS     += -Wall -Woverloaded-virtual -fPIC
DICTCXXFLG   :=
ifdef EXTRAWARN
#FIXME: should be configure'd:
CXXVER       := $(shell g++ --version | head -1 | sed 's/.* \([0-9]\)\..*/\1/')
ifeq ($(CXXVER),4)
CXXFLAGS     += -Wextra -Wno-missing-field-initializers
DICTCXXFLG   := -Wno-strict-aliasing 
endif
endif
LD            = g++
SOFLAGS       = -shared
ifdef I387MATH
CXXFLAGS     += -mfpmath=387
else
CXXFLAGS     += -march=pentium4 -mfpmath=sse
endif
endif

ifeq ($(ARCH),solarisCC5)
# Solaris CC 5.0
CXX           = CC
ifdef DEBUG
  CXXFLAGS    = -g
  LDFLAGS     = -g
  DEFINES     =
else
  CXXFLAGS    = -O
  LDFLAGS     = -O
  DEFINES     = -DNDEBUG
endif
DEFINES      += -DSUNVERS -DHAS_SSTREAM
CXXFLAGS     += -KPIC
LD            = CC
SOFLAGS       = -G
DICTCXXFLG   :=
endif

ifeq ($(CXX),)
$(error $(ARCH) invalid architecture)
endif

ifdef VERBOSE
DEFINES      += -DVERBOSE
endif
ifdef TESTCODE
DEFINES      += -DTESTCODE
endif

CXXFLAGS     += $(DEFINES) $(INCLUDES)
LIBS         += $(ROOTLIBS) $(SYSLIBS)
GLIBS        += $(ROOTGLIBS) $(SYSLIBS)

MAKEDEPEND    = gcc

ifdef PROFILE
CXXFLAGS     += -g -pg
LDFLAGS      += -g -pg
endif

ifndef PKG
PKG           = lib$(CORELIB)
LOGMSG        = "$(PKG) source files"
else
LOGMSG        = "$(PKG) Software Development Kit"
endif
DISTFILE      = $(PKG).tar.gz

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

$(CORELIB):	$(OBJ)
		$(LD) $(LDFLAGS) $(SOFLAGS) -o $@ $^
		@echo "$@ done"

$(MWDCLIB):	$(MOBJ)
		$(LD) $(LDFLAGS) $(SOFLAGS) -o $@ $^ $(CORELIB)
		@echo "$@ done"

$(GEMLIB):	$(GOBJ)
		$(LD) $(LDFLAGS) $(SOFLAGS) -o $@ $^ $(CORELIB)
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
		rm -f $(DISTFILE)
		rm -rf $(PKG)
		mkdir $(PKG)
		cp -p $(SRC) $(HDR) $(LINKDEF) db*.dat Makefile $(PKG)
		cp -p $(MWDCSRC) $(MHDR) $(GEMSRC) $(GHDR) $(PKG)
		gtar czvf $(DISTFILE) --ignore-failed-read \
		 -V $(LOGMSG)" `date -I`" $(PKG)
		rm -rf $(PKG)

.PHONY: all clean realclean srcdist

.SUFFIXES:
.SUFFIXES: .c .cc .cpp .cxx .C .o .d

%.o:	%.cxx
	$(CXX) $(CXXFLAGS) -o $@ -c $<

# FIXME: this only works with gcc
%.d:	%.cxx
	@echo Creating dependencies for $<
	@$(SHELL) -ec '$(MAKEDEPEND) -MM $(INCLUDES) -c $< \
		| sed '\''s%^.*\.o%$*\.o%g'\'' \
		| sed '\''s%\($*\)\.o[ :]*%\1.o $@ : %g'\'' > $@; \
		[ -s $@ ] || rm -f $@'

###

-include $(DEP)
-include $(MDEP)
-include $(GDEP)

