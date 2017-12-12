
# If a special compiler is needed, set PATH and LD_LIBRARY_PATH appropriately
# either in your environment with export or here in a special section for you
SGE_ROOT ?= NOT
ifeq ($(SGE_ROOT), /opt/sge6-2)
  # wigler cluster at CSHL special definitions
  GCC_DIR := /data/software/gcc/4.9.2/rtf
  # GCC_DIR := /data/software/gcc/7.2.0
  # GCC_DIR := /data/software/gcc/6.4.0
  PATH := $(GCC_DIR)/bin:/data/software/local/bin:/bin:/usr/bin
  LD_LIBRARY_PATH := $(GCC_DIR)/lib64
  tcmalloc := defined
  ifdef tcmalloc
    THIRD := /data/software
    tclib := $(THIRD)/gperf/2.1.90/lib
    unwindlib := $(THIRD)/libunwind/0.99-beta/lib
    tclibs := -ltcmalloc
    LD_LIBRARY_PATH := $(LD_LIBRARY_PATH):$(tclib):$(unwindlib)
    LIBFLAGS += -L$(tclib) -L$(unwindlib)
    LIBS += $(tclibs)
  endif
else ifeq ($(SGE_ROOT), /opt/sge)
# NYGC special definitions 
ifeq ($(USER), andrewsp-488)
  # for Peter
  GCC_DIR := /gpfs/commons/home/andrewsp-488/4.9.2
else
  # for Lakshmi
  GCC_DIR := /nfs/sw/gcc/gcc-4.9.2/rtf
endif
  PATH := $(GCC_DIR)/bin:/bin:/usr/bin
  LD_LIBRARY_PATH := $(GCC_DIR)/lib64
else ifeq ($(HOSTNAME), wigstore4.cshl.edu)
  GCC_DIR := /data/software/gcc/4.9.3
  PATH := $(GCC_DIR)/bin:/usr/local/bin:/bin:/usr/bin
  LD_LIBRARY_PATH := $(GCC_DIR)/lib64
endif
