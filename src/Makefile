#======================================================================================/
#**    Author: Joey Azofeifa                                               	       /
#**    Date:   7/8/2015	                                                               /
#**    Description: Makefile for Tfit (Exponentially Modified Gaussian Mixture Model)  /
#======================================================================================/

CXX             = mpic++
CXXFLAGS        = -static-libstdc++ -static-libgcc  -Wno-unused-variable -Wno-non-virtual-dtor -std=c++11 -fopenmp -Wno-write-strings -Wno-literal-suffix -D_GLIBCXX_USE_CXX11_ABI=0 -g
EXEC            = ${PWD}/Tfit
ARCH  			 = getconf LONG_BIT

CPP_FLAGS_32  	= -D32_BIT
CPP_FLAGS_64  	= -D64_BIT

GCCVERSION 	= $(shell ${CXX} -dumpversion)
NU_FIT: main.o load.o split.o model.o across_segments.o template_matching.o \
	read_in_parameters.o model_selection.o error_stdo_logging.o \
	MPI_comm.o  density_profiler.o bootstrap.o prelim_main.o model_main.o select_main.o FDR.o BIC.o ParamWrapper.o old_template_matching.o
	@printf "linking               : "
	@${CXX} ${CXXFLAGS}  ${PWD}/main.o ${PWD}/load.o ${PWD}/model_selection.o \
	${PWD}/split.o ${PWD}/model.o ${PWD}/across_segments.o  \
	${PWD}/template_matching.o ${PWD}/read_in_parameters.o  \
	${PWD}/MPI_comm.o   \
	${PWD}/bootstrap.o ${PWD}/density_profiler.o \
	${PWD}/prelim_main.o ${PWD}/model_main.o ${PWD}/BIC.o ${PWD}/FDR.o  \
	${PWD}/select_main.o  ${PWD}/error_stdo_logging.o ${PWD}/ParamWrapper.o ${PWD}/old_template_matching.o -o ${EXEC} -lmpi
	@cp ${PWD}/Tfit ${PWD}/EMGU
	@printf "done\n"
	@echo "========================================="
	@printf "Successfully compiled! \n\n"

main.o: main.cpp
	@printf "\n=========================================\n"
	@echo "GCC version: "${GCCVERSION}
    
	@printf "main                  : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/main.cpp 
	@printf "done\n"
    
error_stdo_logging.o: error_stdo_logging.cpp
	@printf "error_stdo_logging    : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/error_stdo_logging.cpp 
	@printf "done\n"

BIC.o: BIC.cpp
	@printf "BIC                   : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/BIC.cpp
	@printf "done\n"

select_main.o: select_main.cpp
	@printf "select_main           : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/select_main.cpp 
	@printf "done\n"

model_main.o: model_main.cpp
	@printf "model_main            : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/model_main.cpp 
	@printf "done\n"

prelim_main.o: prelim_main.cpp
	@printf "prelim_main            : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/prelim_main.cpp 
	@printf "done\n"

bootstrap.o: bootstrap.cpp
	@printf "bootstrap             : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/bootstrap.cpp 
	@printf "done\n"

MPI_comm.o: MPI_comm.cpp
	@printf "MPI_comm              : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/MPI_comm.cpp 
	@printf "done\n"

load.o: load.cpp
	@printf "load                  : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/load.cpp 
	@printf "done\n"
    
model_selection.o: model_selection.cpp
	@printf "model_selection       : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/model_selection.cpp 
	@printf "done\n"
    
split.o: split.cpp
	@printf "split                 : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/split.cpp 
	@printf "done\n"

model.o: model.cpp
	@printf "model                 : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/model.cpp 
	@printf "done\n"
    
across_segments.o: across_segments.cpp
	@printf "across_segments       : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/across_segments.cpp 
	@printf "done\n"

template_matching.o: template_matching.cpp
	@printf "template_matching     : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/template_matching.cpp 
	@printf "done\n"

old_template_matching.o: old_template_matching.cpp
	@printf "old_template_matching : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/old_template_matching.cpp 
	@printf "done\n"
    
read_in_parameters.o: read_in_parameters.cpp
	@printf "read_in_parameters    : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/read_in_parameters.cpp 
	@printf "done\n"
    
model_single.o: model_single.cpp
	@printf "model_single          : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/model_single.cpp 
	@printf "done\n"
    
density_profiler.o: density_profiler.cpp
	@printf "density_profiler      : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/density_profiler.cpp 
	@printf "done\n"
    
FDR.o: FDR.cpp
	@printf "FDR                   : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/FDR.cpp
	@printf "done\n"

ParamWrapper.o: ParamWrapper.cpp
	@printf "ParamWrapper          : "
	@${CXX} -c ${CXXFLAGS} ${PWD}/ParamWrapper.cpp 
	@printf "done\n"
	
clean:
	@rm -f *.o
	@rm -f EMGU
	@rm -f Tfit	
	@printf "\nSuccessfully removed binaries!\n\n"
