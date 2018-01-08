
# If a special compiler is needed, set PATH and LD_LIBRARY_PATH appropriately
# either in your environment with export or here in a special section for you
SGE_ROOT ?= NOT
SGE_CLUSTER_NAME ?= NOT
ifeq ($(SGE_CLUSTER_NAME), wigclust)

  # Test code on various compiler versions for warnings and errors
  test_compilers :
	@ unset GCC_DIR && ./make/test_compilers.sh

  # Put online for download
  publish : zip
	scp ~/mumdex.zip mumdex.com:/paa/mumdex.com/
	ssh mumdex.com 'cd /paa/mumdex.com && [ -e mumdex.zip ] && rm -Rf mumdex/ && unzip mumdex.zip'

  ifndef GCC_DIR
    FAST += -march=native -flto
  endif
  # wigler cluster at CSHL special definitions
  GCC_DIR ?= /data/software/gcc/4.9.2/rtf
  # GCC_DIR := /data/software/gcc/4.9.4
  # GCC_DIR := /data/software/gcc/5.5.0
  # GCC_DIR := /data/software/gcc/6.4.0
  # GCC_DIR := /data/software/gcc/7.2.0
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
endif
  PATH := $(GCC_DIR)/bin:$(PATH)
  LD_LIBRARY_PATH := $(GCC_DIR)/lib64:$(LD_LIBRARY_PATH)
endif
