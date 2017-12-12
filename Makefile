#
# Makefile
#
# Makefile for Peter Andrews' source code
#
# Copyright 2017 by Peter Andrews @ CSHL
#

# What to make by default - targets are added below and in module Makefiles
all :

# Shell to use throughout
SHELL := bash

# Which modules to compile and where to find source files
MODULES := core cn convert utility unstable python
VPATH := $(MODULES)
include $(wildcard $(addsuffix /*.mk,$(MODULES)))

# Automatic dependency file creation
DEPDIR := .dep
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.temp_dep
POSTCOMPILE = @ mv -f $(DEPDIR)/$*.temp_dep $(DEPDIR)/$*.dep && touch $@
$(shell mkdir -p $(DEPDIR) > /dev/null)
$(DEPDIR)/%.dep: ;
.PRECIOUS: $(DEPDIR)/%.dep
include $(wildcard $(DEPDIR)/*.dep)

# Compile command and flags
CXX :=  g++
STD := -std=c++11
FAST := -Ofast -flto -m64 -march=native
DEBUG := # -g -ggdb
LDFLAGS	:= $(FAST)
LIBS :=
LIBDIRS := $(addprefix -I ,$(MODULES))

# Excessive compiler warnings are a good thing...
WARN := -Wpedantic -Wall -Wextra -Weffc++ -Wc++11-compat \
	-Wctor-dtor-privacy -Wnarrowing -Wold-style-cast -Woverloaded-virtual \
	-Wsign-promo -Wformat=2 -Wmissing-include-dirs -Wswitch-default \
	-Wswitch-enum -Wunused-parameter -Wuninitialized -Wunknown-pragmas \
	-Wfloat-equal -Wundef -Wshadow -Wlarger-than=10000 \
	-Wframe-larger-than=50000 -Wcast-qual -Wcast-align -Wdate-time \
	-Wenum-compare -Wpacked -Wredundant-decls -Winvalid-pch -Wlong-long \
	-Wvla -Wdisabled-optimization -Wmissing-braces
# Unused warnings to potentially add later:
#
#
# -Wfloat-conversion 
# -Wconversion
# -Wsign-conversion 
# -Wstrict-aliasing=1
# -Wstrict-overflow=5
# -Wsuggest-attribute=pure
# -Wsuggest-attribute=noreturn 
# -Wsuggest-attribute=const 
# -Waggregate-return --- no big deal if so, now in c++14
# -Winline --- no big deal if not inlinable
# -Wpadded --- lots of uncorrectable false positives from lambda structs
# -Wabi --- incompatibilities, not useful
# -Wstack-usage=50000 --- standard headers report unbounded usage
# -Wunsafe-loop-optimizations --- pretty useless?  standard usage triggers it	

# Explicit suffix list : files like these will not be targets in "% :" for speed
.SUFFIXES :
.SUFFIXES : .cpp .h .o .ox .lint .dep

# Do not use shortcut rule - always build .o file (using implicit chaining)
% : %.cpp

# Local machine configuration definitions (edit to fit your setup)
-include make/local.mk

# OS determination and special compilation stuff for each particular OS
UNAME = $(shell sh -c 'uname -s 2> /dev/null || echo NO_UNAME_FOUND')

ifeq ($(UNAME),Linux)
	WARN += -Wnoexcept -Wstrict-null-sentinel -Wdouble-promotion \
		-Wsync-nand -Wtrampolines -Wconditionally-supported \
		-Wlogical-op -Wzero-as-null-pointer-constant \
		-Wvector-operation-performance -Wuseless-cast
	LDFLAGS += -pthread
	ifdef LD_LIBRARY_PATH
		LDFLAGS += -Xlinker -rpath=$(LD_LIBRARY_PATH)
	endif
endif

ifeq ($(UNAME),Darwin)
	WARN += $(COMMON_WARN) -Wno-variadic-macros -Wint-to-void-pointer-cast \
		-Wshorten-64-to-32
	LIBDIRS += -I /opt/X11/include
	LDFLAGS += -L /opt/X11/lib
endif

# For gcc under cygwin
ifeq ($(UNAME),Windows_NT)
#	LIBDIRS += -I /usr/include/
#	LDFLAGS += -L /usr/lib
endif

# Compile
CXXFLAGS =  $(STD) $(FAST) $(DEBUG) $(WARN) $(LIBDIRS)
COMPILE.cpp = $(CXX) $(CXXFLAGS) -c

# Normal compilation and linking pattern rules
%.o : %.cpp
%.o : %.cpp $(DEPDIR)/%.dep
	$(COMPILE.cpp) $(DEPFLAGS) $< -o $@
	$(POSTCOMPILE)
% : %.o ; $(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

# Special rules for X11 compilation and linking
# give .x11o prerequisite if X library needed
%.x11o : %.cpp $(DEPDIR)/%.dep
	$(COMPILE.cpp) $(DEPFLAGS) $< -o $@
	$(POSTCOMPILE)
% : %.x11o ; $(CXX) $(LDFLAGS) $^ -o $@ -lX11 $(LIBS)

# Special rules for gsl compilation and linking
# give .gslo prerequisite if gsl library needed
%.gslo : %.cpp $(DEPDIR)/%.dep
	$(COMPILE.cpp) $(DEPFLAGS) $< -o $@
	$(POSTCOMPILE)
% : %.gslo ; $(CXX) $(LDFLAGS) $^ -o $@ -lgsl -lgslcblas $(LIBS)

# Do not compile gsl programs under Mac OS (remove if gsl available)
ifeq ($(UNAME),Darwin)
%.gslo : %.cpp ; touch $@ && echo skipping compilation of gsl $@ on mac
% : %.gslo ; touch $@ && echo skipping compilation of gsl $@ on mac
endif

# Do not delete .o files when a chain is used
.PRECIOUS : %.o %.gslo %.x11o

# Automatic linting of all source code in all subdirectories
LINTDIR := .lint
LINTCODE := $(wildcard */*.h) $(wildcard */*.cpp)
LINT := ./make/cpplint.py --filter=-readability/streams,-runtime/printf,-build/header_guard,-build/c++11,-runtime/references
LINTGREP := grep -v -e 'Done process' -e 'Total err'
$(shell mkdir -p $(LINTDIR) > /dev/null)
all : lint
lint : $(addprefix $(LINTDIR)/,$(notdir $(LINTCODE:=.lint)))
$(LINTDIR)/%.lint : % ; @ $(LINT) $< 2>&1 | tee >($(LINTGREP) 1>&2) > $@
.PRECIOUS: $(LINTDIR)/%.lint

# Clean all compiled files
clean :
	rm -f *.x11o *.gslo *.o $(PROGRAMS) && rm -Rf $(DEPDIR)/ $(LINTDIR)/
spotless : clean
	rm -f *~ */*~

# Package code as zip in home directory
CODEDIR	= $(notdir $(PWD))
zip :  spotless
	rm -f ~/$(CODEDIR).zip
	(cd .. ; zip -r ~/$(CODEDIR) $(CODEDIR) -x\*.git\* > /dev/null)

# Count code
count : ; @ cat */*.{cpp,h} | sed 's/[ \t]+/ /g' | sort | uniq | wc -l

# Put online for download
publish: zip
	scp ~/mumdex.zip mumdex.com:/paa/mumdex.com/

