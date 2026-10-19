###################################################################################################
#
#  makefile  -  Builds ERMESPETScSolver
#
#  Requires PETSC_DIR and PETSC_ARCH, set in the environment (e.g. in .bashrc):
#     export PETSC_DIR=$HOME/petsc
#     export PETSC_ARCH=arch-complex-M1
#  or given on the command line:
#     make PETSC_DIR=$HOME/petsc PETSC_ARCH=arch-complex-M1
#
#  Compile with:  make    (once for each PETSC_ARCH used)
#
###################################################################################################

# PETSc pkg-config file (named petsc.pc, or PETSc.pc in some versions)
petsc.pc := $(firstword $(wildcard $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/petsc.pc \
                                   $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/PETSc.pc))

ifeq ($(strip $(PETSC_DIR)),)
  $(error PETSC_DIR is not set. Use: export PETSC_DIR=/path/to/petsc  (or make PETSC_DIR=...))
endif
ifeq ($(strip $(PETSC_ARCH)),)
  $(error PETSC_ARCH is not set. Use: export PETSC_ARCH=arch-complex-M1  (or make PETSC_ARCH=...))
endif
ifeq ($(petsc.pc),)
  $(error PETSc not found in $(PETSC_DIR)/$(PETSC_ARCH) (no lib/pkgconfig/petsc.pc). Check PETSC_DIR and PETSC_ARCH)
endif

# Compiler and flags taken from the PETSc installation
CXX      := $(shell pkg-config --variable=cxxcompiler $(petsc.pc))
CXXFLAGS := $(shell pkg-config --variable=cxxflags_extra $(petsc.pc)) $(shell pkg-config --cflags-only-other $(petsc.pc))
CPPFLAGS := $(shell pkg-config --cflags-only-I $(petsc.pc))
LDFLAGS  := $(shell pkg-config --libs-only-L --libs-only-other $(petsc.pc))
LDFLAGS  += $(patsubst -L%, $(shell pkg-config --variable=ldflag_rpath $(petsc.pc))%, $(shell pkg-config --libs-only-L $(petsc.pc)))
LDLIBS   := $(shell pkg-config --libs-only-l $(petsc.pc)) -lm

.PHONY: all clean

all: ERMESPETScSolver

ERMESPETScSolver: ERMESPETScSolver.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f ERMESPETScSolver
