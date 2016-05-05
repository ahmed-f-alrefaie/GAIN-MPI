goal:   GAIN-MPI.x

#
MPICC = mpiicc
FORT = ifort
NVCC = nvcc

NVCCFLAGS = --ptxas-options=-v -O3 -Xptxas -v -lineinfo
#NVCCFLAGS = --ptxas-options=-v -G -g -O0 -Xptxas -v -lineinfo
#Comment on the set of flags for your compiler (GNU or INTEL)
##-------GNU FLAGS-------------------------------------------
#CCFLAGS =   -O3 -fopenmp -march=native
#FFLAGS = -O3 -fopenmp -march=native -ffree-line-length-none
#C_FORT_LIB = -lgfortran


##-------INTEL---------------------------
CCFLAGS =   -O3 -openmp -xHost
FFLAGS = -O3 -openmp -xHost
C_FORT_LIB = -lifcore -limf
#CCFLAGS =   -g -O0 -openmp
#FFLAGS = -g -O0 -openmp
#C_FORT_LIB = -lifcore -limf

#CCFLAGS = -g -O0 -openmp







ifeq ($(PLAT),FERMI)
NVCCFLAGS := $(NVCCFLAGS) -arch sm_21 
else
CCFLAGS := $(CCFLAGS) -DKEPLER
NVCCFLAGS := $(NVCCFLAGS) -arch sm_35 -DKEPLER
PLAT = KEPLER
endif

#C_FORT_LIB = -lifcore -limf


CUDA_LIB = -L$(CUDA_HOME)/lib64/ -lcudart -lcuda -lcublas
LIB     =  $(CUDA_LIB) $(C_FORT_LIB)

CUDA_INC = -I$(CUDA_HOME)/include




###############################################################################

OBJ = BaseProcess.o Util.o Timer.o BaseManager.o TroveDipole.o \
      TroveBasisSet.o TroveInput.o TroveStates.o MultiGpuManager.o \
      EigenVector.o Output.o BasisSet.o Dipole.o GpuManager.o States.o \
      dipole_GPU.o accuracy.o fort_func.o symmetry.o trove_wrappers.o



GAIN-MPI.x:       $(OBJ) main.o
	$(MPICC) -o GAIN-MPI_$(PLAT).x $(OBJ) $(CCFLAGS) main.o $(LIB)  $(CUDA_INC)
#CPP
main.o:       main.cpp
	$(MPICC) -c main.cpp $(CCFLAGS) $(CUDA_INC)

BaseProcess.o: ./common/BaseProcess.cpp
	$(MPICC) -c ./common/BaseProcess.cpp	$(CCFLAGS)

Util.o: ./common/Util.cpp
	$(MPICC) -c ./common/Util.cpp 	$(CCFLAGS)

Timer.o: ./common/Timer.cpp
	$(MPICC) -c ./common/Timer.cpp	$(CCFLAGS)

BaseManager.o: ./common/BaseManager.cpp
	$(MPICC) -c ./common/BaseManager.cpp	$(CCFLAGS)

TroveDipole.o: ./TroveClasses/TroveDipole.cpp
	$(MPICC) -c ./TroveClasses/TroveDipole.cpp	$(CCFLAGS)

TroveBasisSet.o: ./TroveClasses/TroveBasisSet.cpp
	$(MPICC) -c ./TroveClasses/TroveBasisSet.cpp	$(CCFLAGS)

TroveInput.o: ./TroveClasses/TroveInput.cpp
	$(MPICC) -c ./TroveClasses/TroveInput.cpp	$(CCFLAGS)

TroveStates.o: ./TroveClasses/TroveStates.cpp
	$(MPICC) -c ./TroveClasses/TroveStates.cpp	$(CCFLAGS)

MultiGpuManager.o: ./BaseClasses/MultiGpuManager.cpp
	$(MPICC) -c ./BaseClasses/MultiGpuManager.cpp	$(CCFLAGS) $(CUDA_INC)

EigenVector.o: ./BaseClasses/EigenVector.cpp
	$(MPICC) -c ./BaseClasses/EigenVector.cpp	$(CCFLAGS)

Output.o: ./BaseClasses/Output.cpp
	$(MPICC) -c ./BaseClasses/Output.cpp	$(CCFLAGS)

BasisSet.o: ./BaseClasses/BasisSet.cpp
	$(MPICC) -c ./BaseClasses/BasisSet.cpp	$(CCFLAGS)

Dipole.o: ./BaseClasses/Dipole.cpp
	$(MPICC) -c ./BaseClasses/Dipole.cpp	$(CCFLAGS)

GpuManager.o: ./BaseClasses/GpuManager.cpp
	$(MPICC) -c ./BaseClasses/GpuManager.cpp $(CUDA_INC) $(CCFLAGS)

States.o: ./BaseClasses/States.cpp
	$(MPICC) -c ./BaseClasses/States.cpp	$(CCFLAGS)
#FORTRAN
accuracy.o: ./TroveFortran/accuracy.f90
	$(FORT) -c ./TroveFortran/accuracy.f90 $(FFLAGS)
fort_func.o: ./TroveFortran/fort_func.f90 accuracy.o
	$(FORT) -c ./TroveFortran/fort_func.f90 $(FFLAGS)

symmetry.o: ./TroveFortran/symmetry.f90 fort_func.o accuracy.o
	$(FORT) -c ./TroveFortran/symmetry.f90 $(FFLAGS)

trove_wrappers.o: ./TroveFortran/trove_wrappers.f90 symmetry.o fort_func.o accuracy.o
	$(FORT) -c ./TroveFortran/trove_wrappers.f90 $(FFLAGS)
#CUDA
dipole_GPU.o: ./GPU/dipole_GPU.cu
	$(NVCC) -c ./GPU/dipole_GPU.cu $(NVCCFLAGS) $(CUDA_INC)


clean:
	rm $(OBJ) *.mod main.o

