goal:   GAIN-MPI.x

#
MPICC = mpiicc
CCFLAGS =   -O3 -openmp -traceback -xHost

FORT = ifort
FFLAGS = -O3 -openmp -traceback -xHost


NVCC = nvcc
NVCCFLAGS = --ptxas-options=-v -O3 -Xptxas -v -lineinfo

ifeq ($(PLAT),FERMI)
NVCCFLAGS := $(NVCCFLAGS) -arch sm_21 
else
CCFLAGS := $(CCFLAGS) -DKEPLER
NVCCFLAGS := $(NVCCFLAGS) -arch sm_35 -DKEPLER
endif

C_FORT_LIB = -lifcore -limf


CUDA_LIB = -L$(CUDA_HOME)lib64/ -lcudart -lcuda -lcublas
LIB     =  $(CUDA_LIB) $(C_FORT_LIB)

CUDA_INC = -I$(CUDA_HOME)include




###############################################################################

OBJ = BaseProcess.o Util.o Timer.o BaseManager.o TroveDipole.o \
      TroveBasisSet.o TroveInput.o TroveStates.o MultiGpuManager.o \
      EigenVector.o Output.o BasisSet.o Dipole.o GpuManager.o States.o \
      symmetry.o threej.o dipole_GPU.o



GAIN-MPI.x:       $(OBJ) main.o
	$(MPICC) -o GAIN-MPI.x $(OBJ) $(CCFLAGS) main.o $(LIB)  $(CUDA_INC)
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
symmetry.o: ./common/symmetry.f90 threej.o
	$(FORT) -c ./common/symmetry.f90 $(FFLAGS)

threej.o: ./common/threej.f90
	$(FORT) -c ./common/threej.f90 $(FFLAGS)
#CUDA
dipole_GPU.o: ./GPU/dipole_GPU.cu
	$(NVCC) -c ./GPU/dipole_GPU.cu $(NVCCFLAGS) $(CUDA_INC)


clean:
	rm $(OBJ) *.mod main.o

