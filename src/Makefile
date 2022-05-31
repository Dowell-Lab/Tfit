# Makefile for Tfit, rewritten July 2021 

VERSION = 1.0
VDATE = JUL 2021

# Compiler choices
#CXX             = /usr/lib64/mpich/bin/mpicxx		# mpich
CXX             = mpicxx				# openmpi or mpich
CXXFLAGS        = -O2 -static-libstdc++ -static-libgcc -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11 -fopenmp -Wno-write-strings -Wno-literal-suffix -D_LGIBCXX_USE_CXX1_ABI=0 -g -pg

EXEC            = ${PWD}/Tfit

ARCH  		= getconf LONG_BIT
CPP_FLAGS_32  	= -D32_BIT
CPP_FLAGS_64  	= -D64_BIT

GCCVERSION 	= $(shell ${CXX} -dumpversion)

PROGS = Tfit argTest
OBJ = load.o split.o model.o across_segments.o template_matching.o \
      read_in_parameters.o model_selection.o error_stdo_logging.o\
      MPI_comm.o density_profiler.o bootstrap.o bidir_main.o model_main.o \
      select_main.o FDR.o BIC.o EMGParams.o helper.o Bedfile.o \
			Intervals.o Regions.o Data.o ITree.o Distro.o Models.o ModelSets.o \
			EMalg.o ModelSupport.o
SRC = $(OBJ:.o=.cpp)

OBJ2 = parse_args.o
SRC2 = $(OBJ2:.o=.cpp)
CFLAGS = -g -std=c++11

### Build instructions

Tfit: main.o ${OBJ}
	@printf "linking           : "
	@${CXX} -o ${EXEC} ${CXXFLAGS} main.o ${OBJ} -lmpi
	@printf "done\n"
	@echo "========================================="
	@printf "Tfit version: "${VERSION}
	@printf " successfully compiled \n\n"

argTest:
	echo "Temporary Argument Stuff"
	g++ -c ${CFLAGS} ${PWD}/parse_args.cpp 
	g++ -c ${CFLAGS} ${PWD}/arg_main.cpp 
	g++ -o argTest ${CFLAGS} arg_main.o parse_args.o 

testmain: MainTesting.o ${OBJ}
	echo "Temporary Testing"
	@${CXX} -o Testmain ${CXXFLAGS} MainTesting.o ${OBJ} -lmpi

%: 
	@${CXX} -c ${CXXFLAGS} ${PWD}/$*.cpp

clean:
	@rm -f *.o $(PROGS)
	@printf "\nsuccessfully removed binaries\n\n"
