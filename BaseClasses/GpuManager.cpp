#include "GpuManager.h"


GpuManager(int pgpu_id,int pmax_degen,int nprocs) : BaseProcess(), BaseManager(), gpu_id(pgpu_id), max_degen(pmax_degen),Nprocs(nprocs){



}
void AllocateVectors(int nsizemax,int dimenmax){



}
void InitializeAndTransferConstants(int jmax,int ){
	//set the device
	cudaSetDevice(gpu_id);
	cudaFree(0);
	cudaCheckError("Wake Device");




}

void TransferBasisSet(const BasisSet* basisSet);
void TransferInflation(int* icontr_, int* ijterm);
void TransferDipole(Dipole* dipole_);

	

void ExecuteHalfLs(double* half_ls,double * vector, int N,int indI,int indF,int idegI,int jF,int igammaI);
void ExecuteDotProduct(double* int indI,int indF,int idegI,int jF,int igammaI);
