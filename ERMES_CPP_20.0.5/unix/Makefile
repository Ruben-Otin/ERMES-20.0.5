#==============================================================================#
#  ERMES 20.0 - Linux Makefile                                                 #
#                                                                              #
#  Hand-written replacement for the cbp2make-generated Makefile, which is     #
#  kept as Makefile_codeblocks.                                                #
#                                                                              #
#  Usage:                                                                      #
#    make -j$(nproc)              optimized build ( -O3 ) -> ./ERMES           #
#    make -j NATIVE=1             also tune for the build machine CPU          #
#                                 ( only if you run on the same CPU type )    #
#    make -j STATIC=1             static executable, portable to older Linux  #
#    make -j METIS=0              build without the bundled METIS             #
#    make -j DEBUG=1              debug build ( -O0 -g, bounds-checked STL )  #
#    make -j SANITIZE=1           debug build with AddressSanitizer / UBSan    #
#    make clean                                                                #
#                                                                              #
#  Changes with respect to the previous Makefile:                              #
#   - Warnings are enabled for ERMES sources ( previously -w hid all of them ) #
#     and -fpermissive is no longer needed.                                   #
#   - Explicit -std=c++11: GCC >= 11 defaults to C++17, which the old code    #
#     did not compile with.                                                    #
#   - Header dependencies are tracked ( -MMD ): editing a .h rebuilds the     #
#     files that include it ( the old Makefile did not rebuild them ).        #
#   - Objects go to build/ ( no spaces in paths ).                             #
#==============================================================================#

CXX      ?= g++
CC       ?= gcc

EXE      := ERMES
BUILDDIR := build

ROOT     := ..

#------------------------------------------------------------------------------
# Flags
#------------------------------------------------------------------------------
STD      := -std=c++11

# Boost is no longer needed ( std::shared_ptr is used instead of boost::shared_ptr )
INC      := -I$(ROOT)/external_libraries/ \
            -I$(ROOT)/external_libraries/gidpost/ \
            -I$(ROOT)/external_libraries/gidpost/win/ \
            -I$(ROOT)/external_libraries/flex++bison++/ \
            -I$(ROOT)/external_libraries/flex++bison++/flex_bison/

# Warnings for ERMES code. The silenced ones are pervasive in the legacy code
# base ( virtual interfaces with unused arguments, signed/unsigned loops, member
# init order ) and do not indicate bugs.
WARN     := -Wall -Wextra -Wno-unused-parameter -Wno-sign-compare -Wno-reorder \
            -Wno-overloaded-virtual -Wno-deprecated-copy -Wno-deprecated-declarations \
            -Wno-unused-variable -Wno-unused-but-set-variable -Wno-shadow -Wno-misleading-indentation -Wno-write-strings

ifeq ($(DEBUG),1)
  OPT    := -O0 -g -D_GLIBCXX_ASSERTIONS
else ifeq ($(SANITIZE),1)
  OPT    := -O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined
else
  OPT    := -O3 -DNDEBUG
  ifeq ($(NATIVE),1)
    OPT  += -march=native
  endif
endif

OMP      := -fopenmp

#------------------------------------------------------------------------------
# Optional libraries for the internal direct solver ( open source )
#
#   make BLAS=openblas            optimized dense kernels with OpenBLAS ( -lopenblas )
#   make BLAS=blis                optimized dense kernels with BLIS     ( -lblis )
#   make BLAS=other BLAS_LIB="-L/path/lib -lmyblas"   any BLAS with the CBLAS interface
#   ( METIS is bundled and included by default, see below )
#
#   With BLAS=openblas or BLAS=blis ERMES sets the BLAS to one thread itself ( it distributes
#   the work over its own OpenMP threads ). With BLAS=other set e.g. MKL_NUM_THREADS=1 yourself.
#------------------------------------------------------------------------------
DEFS    :=
EXTLIBS :=

ifneq ($(BLAS),)
  DEFS += -DERMES_USE_BLAS
  ifeq ($(BLAS),openblas)
    DEFS += -DERMES_USE_OPENBLAS
    BLAS_LIB ?= -lopenblas
  endif
  ifeq ($(BLAS),blis)
    DEFS += -DERMES_USE_BLIS
    BLAS_LIB ?= -lblis
  endif
  EXTLIBS += $(BLAS_LIB)
endif

# METIS ( nested dissection ordering for the direct solver ) is included by default: its sources are
# bundled in external_libraries/metis and external_libraries/GKlib and compiled together with ERMES.
#   make METIS=0                       build without METIS
#   make METIS_DIR=/path [METIS_LIB=]  use an installed METIS instead of the bundled sources
METIS ?= 1

METIS_BUNDLED := $(wildcard $(ROOT)/external_libraries/metis/libmetis/*.c)

METIS_C_SRC :=
ifeq ($(METIS),1)
  DEFS += -DERMES_USE_METIS
  ifneq ($(METIS_DIR),)
    METIS_LIB ?= -lmetis
    INC       += -I$(METIS_DIR)/include
    EXTLIBS   += -L$(METIS_DIR)/lib $(METIS_LIB)
  else
    ifeq ($(METIS_BUNDLED),)
      $(error METIS sources not found in external_libraries/metis: use METIS=0 or METIS_DIR=/path)
    endif
    INC         += -I$(ROOT)/external_libraries/metis/include
    METIS_C_SRC := $(METIS_BUNDLED) $(wildcard $(ROOT)/external_libraries/GKlib/src/*.c)
    EXTLIBS     += -lm
  endif
endif

CXXFLAGS_ERMES := $(STD) $(OPT) $(OMP) $(WARN) $(DEFS) $(INC) -MMD -MP
CXXFLAGS_EXT   := $(STD) $(OPT) $(OMP) -w   $(INC) -MMD -MP
CFLAGS_EXT     :=        $(OPT)        -w   $(INC) -MMD -MP

LDFLAGS  += $(OMP)

# Fully static executable ( glibc, libstdc++, libgomp linked in ): runs on any x86-64 Linux,
# including clusters with an older system library than the build machine.
ifeq ($(STATIC),1)
  LDFLAGS += -static
endif
ifeq ($(SANITIZE),1)
  LDFLAGS += -fsanitize=address,undefined
endif

#------------------------------------------------------------------------------
# Sources ( sources/update.cpp is not part of the build: it is legacy code
# that does not compile )
#------------------------------------------------------------------------------
ERMES_SRC := $(wildcard $(ROOT)/ERMES/*.cpp) \
             $(wildcard $(ROOT)/ERMES/Utils/*.cpp) \
             $(wildcard $(ROOT)/linear_solvers/*.cpp) \
             $(filter-out $(ROOT)/sources/update.cpp, $(wildcard $(ROOT)/sources/*.cpp))

EXT_CXX_SRC := $(ROOT)/external_libraries/gidpost/source/gidpost.cpp \
               $(ROOT)/external_libraries/gidpost/source/gidpostInt.cpp

EXT_C_SRC   := $(ROOT)/external_libraries/gidpost/source/zlibint1.c \
               $(ROOT)/external_libraries/gidpost/source/zlibint2.c

obj = $(patsubst $(ROOT)/%,$(BUILDDIR)/%,$(1))

ERMES_OBJ   := $(call obj,$(ERMES_SRC:.cpp=.o))
EXT_CXX_OBJ := $(call obj,$(EXT_CXX_SRC:.cpp=.o))
EXT_C_OBJ   := $(call obj,$(EXT_C_SRC:.c=.o))

METIS_OBJ   := $(call obj,$(METIS_C_SRC:.c=.o))

ALL_OBJ     := $(ERMES_OBJ) $(EXT_CXX_OBJ) $(EXT_C_OBJ) $(METIS_OBJ)

#------------------------------------------------------------------------------
# Rules
#------------------------------------------------------------------------------
.PHONY: all clean

all: $(EXE)

$(EXE): $(ALL_OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS) $(EXTLIBS)

$(ERMES_OBJ): $(BUILDDIR)/%.o: $(ROOT)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS_ERMES) -c $< -o $@

$(EXT_CXX_OBJ): $(BUILDDIR)/%.o: $(ROOT)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS_EXT) -c $< -o $@

$(EXT_C_OBJ): $(BUILDDIR)/%.o: $(ROOT)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS_EXT) -c $< -o $@

# Bundled METIS / GKlib ( C99, same settings as METIS's own Linux build, 32-bit indices )
METIS_CFLAGS := $(OPT) -std=c99 -fno-strict-aliasing -w -DLINUX -D_POSIX_C_SOURCE=200809L -D_FILE_OFFSET_BITS=64 \
                -DNDEBUG -DNDEBUG2 -DHAVE_EXECINFO_H -DHAVE_REGEX_H -DHAVE_GETLINE \
                -I$(ROOT)/external_libraries/metis/include -I$(ROOT)/external_libraries/GKlib/include \
                -I$(ROOT)/external_libraries/metis/libmetis -MMD -MP

$(METIS_OBJ): $(BUILDDIR)/%.o: $(ROOT)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(METIS_CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILDDIR) $(EXE)

-include $(ALL_OBJ:.o=.d)

#------------------------------------------------------------------------------
# Unit tests ( no input files needed )
#------------------------------------------------------------------------------
.PHONY: test

$(BUILDDIR)/test_matrix_operations: $(ROOT)/tests/test_matrix_operations.cpp $(BUILDDIR)/linear_solvers/complex_matrix_operations.o
	@mkdir -p $(dir $@)
	$(CXX) $(STD) -O2 $(OMP) $(INC) $^ -o $@

$(BUILDDIR)/test_direct_solver: $(ROOT)/tests/test_direct_solver.cpp $(BUILDDIR)/linear_solvers/direct_complex_solver.o \
                                 $(BUILDDIR)/linear_solvers/multifrontal_complex_solver.o \
                                 $(BUILDDIR)/linear_solvers/complex_solver.o $(BUILDDIR)/linear_solvers/complex_matrix_operations.o \
                                 $(BUILDDIR)/ERMES/Utils/MatrixFileWritting.o $(METIS_OBJ)
	@mkdir -p $(dir $@)
	$(CXX) $(STD) -O2 $(OMP) $(DEFS) $(INC) $^ -o $@ $(EXTLIBS)

test: $(BUILDDIR)/test_matrix_operations $(BUILDDIR)/test_direct_solver
	OMP_DYNAMIC=false ./$(BUILDDIR)/test_matrix_operations
	./$(BUILDDIR)/test_direct_solver
