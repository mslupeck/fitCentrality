ALIROOT_SO_DIR = /home/ubuntu1804/mss/software/alice/sw/ubuntu1804_x86-64/AliRoot/latest/lib
#---------------=============== ROOT  RELATED ===============---------------
ROOTCONFIG   := root-config

ARCH         := $(shell $(ROOTCONFIG) --arch)
PLATFORM     := $(shell $(ROOTCONFIG) --platform)
ALTCC        := $(shell $(ROOTCONFIG) --cc)
ALTCXX       := $(shell $(ROOTCONFIG) --cxx)
ALTF77       := $(shell $(ROOTCONFIG) --f77)
ALTLD        := $(shell $(ROOTCONFIG) --ld)

#CXX           =
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o # keep whitespace after "-o"

ifeq (debug,$(findstring debug,$(ROOTBUILD)))
OPT           = -g
OPT2          = -g
else
ifneq ($(findstring debug, $(strip $(shell $(ROOTCONFIG) --config))),)
OPT           = -g
OPT2          = -g
else
OPT           = -O
OPT2          = -O2
endif
endif

ROOTCFLAGS   := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS  := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS     := $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS    := $(shell $(ROOTCONFIG) --glibs)
HASTHREAD    := $(shell $(ROOTCONFIG) --has-thread)
ROOTDICTTYPE := $(shell $(ROOTCONFIG) --dicttype)
#NOSTUBS      := $(shell $(ROOTCONFIG) --nostubs)
ROOTCINT     := rootcint

# Stub Functions Generation
ifeq ($(NOSTUBS),yes)
   ROOTCINT = export CXXFLAGS="$(CXXFLAGS)"; $(ROOTSYS)/core/utils/src/rootcint_nostubs.sh -$(ROOTDICTTYPE)
endif

ifeq ($(ARCH),linuxx8664gcc)
# Linux with egcs, gcc 2.9x, gcc 3.x
CXX           = g++
CXXFLAGS      = $(OPT2) -Wall -fPIC
LD            = g++
#LDFLAGS       = $(OPT2) -L/usr/local/lib -lnsl -lz -lm -L$(ALIROOT_SO_DIR) -lALIGN -lANALYSISaliceBase -lANALYSISalice -lANALYSIScalib -lANALYSIS -lSTAT -lSTEERBase -lSTEER -lAOD -lBCM -lCDB -lESD -lFITbase -FITcalib -FITrec -FITsim
LDFLAGS       = $(OPT2) -L/usr/local/lib -lnsl -lz -lm
SOFLAGS       = -shared
endif

CXXFLAGS     += $(ROOTCFLAGS)
LDFLAGS      += $(ROOTLDFLAGS)
LIBS          = $(ROOTLIBS) $(SYSLIBS)
GLIBS         = $(ROOTGLIBS) $(SYSLIBS)

#GLIBS        += -L$(ALIROOT_SO_DIR) -lSTEER -lSTEERBase -lTHijing -lEG
GLIBS        += -L$(ALIROOT_SO_DIR) -lSTEER -lSTEERBase -lTHijing -lEG -lAliPythia6 -lAliPythia8 -lpythia8 -lpythia6

ifneq ($(ALTCC),)
   CC  = $(ALTCC)
endif
ifneq ($(ALTCXX),)
   CXX = $(ALTCXX)
endif
ifneq ($(ALTF77),)
   F77 = $(ALTF77)
endif
ifneq ($(ALTLD),)
   LD  = $(ALTLD)
endif

ifneq ($(findstring g++, $(CXX)),)
GCC_MAJOR := $(shell $(CXX) -dumpversion 2>&1 | cut -d'.' -f1)
GCC_MINOR := $(shell $(CXX) -dumpversion 2>&1 | cut -d'.' -f2)
endif

CFLAGS	= -Wall -ggdb $(CXXFLAGS) $(LDFLAGS) $(ROOTCFLAGS) $(ROOTLDFLAGS) $(ROOTLIBS) $(ROOTGLIBS) -I/home/mss/software/alice/sw/ubuntu1804_x86-64/AliRoot/latest/include
#---------------=============== ROOT  RELATED ===============---------------

#---------------=============== MAIN MAKEFILE ===============---------------
CC			= g++
TARGET		= fitCentrality 
SOTARGET	= libTreeStorage.so
ROOTDICT	= fitCentrality__Dict.cpp
ROOTDICTSRC	= TreeStorage.h LinkDef.h
SOURCE		= $(ROOTDICT) TreeStorage.cpp fitCentrality.cpp

INCLUDE	= -I./ $(AN_INCLUDE_DIR)  
## end more includes

OBJ=$(join $(addsuffix obj/, $(dir $(SOURCE))), $(notdir $(SOURCE:.cpp=.o))) 

## Fix dependency destination to be ../.dep relative to the src dir
DEPENDS=$(join $(addsuffix .dep/, $(dir $(SOURCE))), $(notdir $(SOURCE:.cpp=.d)))

## Default rule executed
all: $(SOTARGET) $(TARGET)
	@true

## Clean Rule
clean:
	@echo "Removing .o files"
	@-rm -f $(TARGET) $(SOTARGET) $(OBJ) $(DEPENDS)
	@echo "Removing root dictionary files"
	@-rm -f $(ROOTDICT) $(ROOTDICT:.cpp=.h)
	@-rm -f AutoDict_*cxx*
	@echo "Removing obj folders"
	@-rmdir obj
	@echo "Removing .dep folders"
	@-rmdir .dep

## Rule for making the .so library
$(SOTARGET): $(OBJ)
	@echo "============="
	@echo "Linking the .so lib $@"
	@echo "============="
	@$(CC) $(SOFLAGS) -o $@ $^ $(GLIBS)
	@echo -- Lib finished --

## Rule for making the actual target
$(TARGET): $(OBJ)
	@echo "============="
	@echo "Linking the target $@"
	@echo "============="
	@$(CC) $(CFLAGS) -o $@ $^ $(GLIBS)
	@echo -- Link finished --

## Generic compilation rule
%.o : %.cpp
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -c $< -o $@


## Rules for object files from .cpp files
## Object file for each file is put in obj directory
obj/%.o : %.cpp
	@mkdir -p $(dir $@)
	@echo "============="
	@echo "Compiling $<"
	@$(CC) $(CFLAGS) -c $< -o $@

## Make dependancy rules
.dep/%.d: %.cpp
	@mkdir -p $(dir $@)
	@echo "============="
	@echo Building dependencies file for $*.o
	@$(SHELL) -ec '$(CC) -M $(CFLAGS) $< | sed "s^$*.o^obj/$*.o^" > $@'

## Make root dictionary
$(ROOTDICT): $(ROOTDICTSRC)
	@echo "============="
	@echo "Generating dictionary $<"
	@$(ROOTCINT) -f $(ROOTDICT) -c $(ROOTDICTSRC)

## Include the dependency files
-include $(DEPENDS)
#---------------=============== MAIN MAKEFILE ===============---------------
