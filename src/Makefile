# Makefile for Tfit, rewritten July 2021 

VERSION = 1.0
VDATE = JUL 2021

# Compiler choices
CXX             = mpic++
CXXFLAGS        = -static-libstdc++ -static-libgcc -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11 -fopenmp -Wno-write-strings -Wno-literal-suffix -D_LGIBCXX_USE_CXX1_ABI=0 -g

EXEC            = ${PWD}/Tfit

ARCH  		= getconf LONG_BIT
CPP_FLAGS_32  	= -D32_BIT
CPP_FLAGS_64  	= -D64_BIT

GCCVERSION 	= $(shell ${CXX} -dumpversion)

PROGS = Tfit
OBJ = load.o split.o model.o across_segments.o template_matching.o \
      read_in_parameters.o model_selection.o error_stdo_logging.o\
      MPI_comm.o density_profiler.o bootstrap.o bidir_main.o model_main.o \
      select_main.o FDR.o BIC.o
SRC = $(OBJ:.o=.cpp)

### Build instructions

Tfit: main.o ${OBJ}
	@printf "linking           : "
	@${CXX} -o ${EXEC} ${CXXFLAGS} main.o ${OBJ} -lmpi
	@printf "done\n"
	@echo "========================================="
	@printf "Tfit version: "${VERSION}
	@printf " successfully compiled \n\n"

%: 
	@${CXX} -c ${CXXFLAGS} ${PWD}/$*.cpp

clean:
	@rm -f *.o $(PROGS)
	@printf "\nsuccessfully removed binaries\n\n"
