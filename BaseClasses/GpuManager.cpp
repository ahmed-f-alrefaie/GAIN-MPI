#include "GpuManager.h"


GpuManager::GpuManager(int pgpu_id,int nprocs) : BaseProcess(), BaseManager(), gpu_id(pgpu_id), Nprocs(nprocs){



}
void GpuManager::AllocateVectors(int nsizemax,int dimenmax){



}
void GpuManager::InitializeAndTransferConstants(int jmax,int sym_repres,int pmax_degen){
	//set the device

	Log("Initializing GPU device...",gpu_id);

	cudaSetDevice(gpu_id);
	cudaFree(0);
	cudaCheckError("Wake Device");
	

	//Lets get some data goin:
	cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, gpu_id);
	
	//Lets get some gpu_info.
	InitializeMemory(size_t(double(devProp.totalGlobalMem)*0.95));
	Log("Gpu ID %d : Total Memory: %12.6f\n",gpu_id,
	


	//Some constants we might need
	MaxDegen = pmax_degen;
	SymNrepres = sym_nrepres;
	Jmax = jmax;

	//Copy constants to the GPU
	cudaMemcpyToSymbol (g_maxdegen,&MaxDegen, sizeof(int) );		
	cudaMemcpyToSymbol (c_jmax,&Jmax , sizeof(int) );		
	cudaMemcpyToSymbol ( g_symnrepres,&SymNrepres, sizeof(int) );	
	CheckCudaError("Copy Constants");


	//Allocate streams
	//Half_linestrength
	for(int i = 0; i < (2*Jmax+1); i++){
		half_ls_stream.push_back(NULL);
		cudaStreamCreate(&half_ls_stream.back());
	}
	

	//DotProduct
	for(int i =0; i < Nprocs; i++){
		dot_product_omp_stream.push_back(NULL);
		cudaStreamCreate(&dot_product_omp_stream.back());
	}


	Log("..Done!\n");


}

void GpuManager::TransferBasisSet(const BasisSet* basisSet){
	cudaSetDevice(gpu_id);

}

void AllocateVectors(int nsizemax,int dimenmax){
	cudaSetDevice(gpu_id);
	Log("Allocating Vectors in GPU %d...",gpu_id);
	
}


void GpuManager::TransferInflation(int* icontr_, int* ijterm,int dimen,std::vector<int> Ntot,std::vector<int> sym_degen){
	cudaSetDevice(gpu_id);
}
void GpuManager::TransferDipole(Dipole* dipole_,int which_block){
	cudaSetDevice(gpu_id);
}

	

void GpuManager::ExecuteHalfLs(double* half_ls,double * vector, int N,int indI,int indF,int idegI,int jF,int igammaI){
	cudaSetDevice(gpu_id);

}
void GpuManager::ExecuteDotProduct(double* int indI,int indF,int idegI,int jF,int igammaI){
	cudaSetDevice(gpu_id);
}
