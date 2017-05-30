#################################################
# Makefile for FIRST
# February 17, 2005. BMH
#
#
#
# 2007.06.05 - BMH - Added script to do full
# dependency checking (depend.sh). Modified
# Makefile to use "modules" for compartmentalized
# source code, and recursive makefiles for
# independent ancilliary programs.
#
# 2007.08.01 - BMH - Minor changes to reflect new
# file names. Removed rule for depends file, and
# switched to just use the *.d files.
#
# 2006.02.11 - SAM - Cleaned up a bit of redundancy
# and added a `depends' target so that modifications
# to a header will no longer force a recompile of
# every source file (now, only those files that depend
# on the header should be rebuilt)
#
# 2005.11.08 - BMH - Modified file to use Main.cpp
# instead of FIRST5.cpp. This reflects the new
# CVS development trunk named FIRST (major and minor
# releases are now indexed by CVS revision tags).
##################################################

FIRST := FIRST
DIAGSTD := diagstd
CUR_VER = FIRST-dev

# List external programs that need to be compiled
# independantly.
EXT_PROGS := hbdilute

# List the module directories.
MODULES := Assembly TIMME froda2 XML

INC := $(patsubst %, -I%/, $(MODULES))
INC += -I../include/

LIBS :=
SRC := $(wildcard *.cpp)

FORT = gfortran
CPP = g++
CPP_OPTIMIZATION = -O3
#DEBUG = -g -pg
EXTRA = -Wall
CPPFLAGS = $(INC) ${CPP_OPTIMIZATION} $(DEBUG) ${EXTRA} -ansi -fPIC

INSTALL_HOME = $(CURDIR)/..
DIAGSTD_HOME = $(CURDIR)/../..

FIRST_ROOT = -DFIRST_ROOT='"${INSTALL_HOME}"'

# Include the makefile info from each of the modules. These
# need to be named "module.mk".
include $(patsubst %, %/module.mk, $(MODULES))

# Get the full list of object files to be linked into
# FIRST from all the modules.
OBJLIST := $(patsubst %.cpp, %.o, $(filter %.cpp, $(SRC)))
DEPLIST := $(patsubst %.cpp, %.d, $(filter %.cpp, $(SRC)))
#EXCLUDED = hbdilute.o
#OBJLIST := $(filter-out $(EXCLUDED), $(OBJLIST))

# The default make will depend on FIRST and any other
# INDEPENDENT executables.
all: FIRST external_progs DIAGSTD

# Compile FIRST
FIRST: $(OBJLIST)
	$(CPP) -o $@ $(OBJLIST) $(LIBS) $(CPPFLAGS)

DIAGSTD: diagstd.f
	$(FORT) -o diagstd diagstd.f

%.o: %.cpp
	${CPP} ${CPPFLAGS} ${FIRST_ROOT} -c $< -o $@

-include $(OBJLIST:.o=.d)

%.d: %.cpp
	./depend.sh ${CPP} `dirname $<` $(CPPFLAGS) $< > $@

%.d: %.h
	./depend.sh ${CPP} `dirname $<` $(CPPFLAGS) $< > $@

clean:
	@for i in $(OBJLIST); do \
	\rm -f $$i; done
	@for i in $(DEPLIST); do \
	\rm -f $$i; done
	@for i in $(EXT_PROGS); do \
        (cd ../$$i; $(MAKE) clean); done

############################################################
# Any independent executable programs (there is only
# one at this time, and hopefully it will be gone soon).
# Should have thier own makefile and be compiled
# recursively.
external_progs:
	@for i in $(EXT_PROGS); do \
        echo "Compiling in $$i..."; \
        (cd ../$$i; $(MAKE) $(MFLAGS)); done

print:
	@echo $(DEPLIST)
