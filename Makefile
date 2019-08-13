#
# Makefile
#
# Makefile for Peter Andrews' source code
#
# Copyright 2019 by Peter Andrews @ CSHL
#

# What to make by default - targets are added below and in module Makefiles
all :

# Shell to use throughout
# SHELL = $(warning [$@])bash
SHELL = bash

# Which modules to compile and where to find source files
include project.mk
VPATH := $(MODULES) $(INCLUDES)
include $(wildcard $(addsuffix /*.mk,$(MODULES)))

# Automatic dependency file creation and inclusion
DEPDIR := .dep
DEPFLAGS = -MT $@ -MMD -MP -MF $(DEPDIR)/$*.temp_dep
POSTCOMPILE = @ mv -f $(DEPDIR)/$*.temp_dep $(DEPDIR)/$*.dep && touch $@
$(shell mkdir -p $(DEPDIR) > /dev/null)
$(DEPDIR)/%.dep: ;
.PRECIOUS: $(DEPDIR)/%.dep
include $(wildcard $(DEPDIR)/*.dep)

# Compilation variables
STD := -std=c++11

# debug = 1
ifdef debug
  FAST :=  -g -ggdb
else
  FAST := -Ofast -march=native -DNDEBUG
  # For gcc under cygwin
  ifneq (,$(findstring /cygdrive/,$(PATH)))
    FAST += -Wa,-mbig-obj
  else
    LTO := -flto
    FAST += $(LTO)
  endif
endif
WARN := -Wall
LIBS := -pthread
INCFLAGS := $(addprefix -I ,$(INCLUDES))
LDFLAGS := $(LTO)

# OS determination and special compilation stuff for each particular OS
UNAME := $(shell sh -c 'uname -s 2> /dev/null || echo NO_UNAME_FOUND')
ifeq ($(UNAME),Darwin)
  CXX := clang++
  INCFLAGS += -I /opt/X11/include
  LDFLAGS += -L /opt/X11/lib
else
  CXX ?= g++
endif

# Local machine configuration definitions (possibly edit to fit your setup)
-include make/local.mk

# Build library paths into the executables, if needed
ifeq ($(CXX),g++)
  ifdef LD_LIBRARY_PATH
    LDFLAGS += -Xlinker -rpath=$(LD_LIBRARY_PATH)
  endif
endif

# Compile
CXXFLAGS = $(STD) $(FAST) $(WARN) $(INCFLAGS)
COMPILE.cpp = $(CXX) $(CXXFLAGS) -c

# Explicit suffix list : files like these will not be targets in "% :" for speed
.SUFFIXES :
#.SUFFIXES : .cpp .h .o .x11o .gslo .eo .lint .dep
.SUFFIXES : .cpp .h .o .ox .lint .dep

# Do not use shortcut rule - always build .o file (using implicit chaining)
% : %.cpp

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

# Special rules for eigen compilation and linking
# give .eo prerequisite if eigen headers are needed
%.eo : %.cpp $(DEPDIR)/%.dep
	$(COMPILE.cpp) $(DEPFLAGS) $< -o $@
	$(POSTCOMPILE)
% : %.eo ; $(CXX) $(LDFLAGS) $^ -o $@ $(LIBS)

# Do not delete .o files when a chain is used
.PRECIOUS : %.o %.gslo %.x11o %.eo

# Automatic linting of all source code in all subdirectories
LINTDIR := .lint
LINTCODE := $(wildcard $(addsuffix /*.h,$(MODULES))) $(wildcard $(addsuffix /*.cpp,$(MODULES)))
LINT := ./make/cpplint.py --filter=-readability/streams,-runtime/printf,-build/header_guard,-build/c++11,-runtime/references
LINTGREP := grep -v -e 'Done process' -e 'Total err'
$(shell mkdir -p $(LINTDIR) > /dev/null)
all : lint
lint : $(addprefix $(LINTDIR)/,$(notdir $(LINTCODE:=.lint)))
$(LINTDIR)/%.lint : % ; @ $(LINT) $< 2>&1 | tee >($(LINTGREP) 1>&2) > $@
.PRECIOUS: $(LINTDIR)/%.lint

# Clean all compiled files
clean :
	rm -f *.x11o *.gslo *.eo *.o $(PROGRAMS)
	rm -Rf $(DEPDIR)/ $(LINTDIR)/ python/build python/dist
spotless : clean
	rm -f *~ */*~

# Package code as zip in home directory
CODEDIR = $(notdir $(PWD))
zip :  spotless
	rm -f ~/$(CODEDIR).zip
	(cd .. ; zip -r ~/$(CODEDIR) $(CODEDIR) -x\*.git\* > /dev/null)

# Count code
count : ; @ cat */*.{cpp,h} | sed 's/[ \t]+/ /g' | sort | uniq | wc -l

# setup binary directory
bin : all
	@mkdir -p ~/bin
	@make fbin
fbin :
	@rsync $(shell find $(shell pwd)/ -perm -u=x -type f | \
		grep -v -e \.git -e /python/ -e '\'~$\'') ~/bin/
	@echo bin transfer done
