#include "GpuManager.h"


GpuManager::GpuManager(int pgpu_id,int nprocs) : BaseProcess(), BaseManager(), gpu_id(pgpu_id), Nprocs(nprocs){



}
void GpuManager::AllocateVectors(int nsizemax,int dimenmax){



}

void GpuManager::AllocateGpuMemory(void** mem, size_t size){
		if(cudaSuccess != cudaMalloc(mem,size)){
			CheckCudaError("Memory Allocation");
		}
		TrackMemory(size);
}

void GpuManager::AllocatePinnedMemory(void** mem, size_t size){
		if(cudaSuccess != cudaMallocHost(mem,size){
			CheckCudaError("Memory Allocation");
		}
		BaseManager::TrackGlobalMemory(size);
}

void GpuManager::TransferToGpu(void* dst,void* src,size_t size){
		if(cudaSuccess != cudaMemcpy(dst,src,size,cudaMemcpyHostToDevice)){
			CheckCudaError("Memory Transfer H -> D");
		}
}

void GpuManager::TransferToHost(void* dst,void* src,size_t size){
		if(cudaSuccess != cudaMemcpy(dst,src,size,cudaMemcpyDeviceToHost)){
			CheckCudaError("Memory Transfer H -> D");
		}
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
	Log("Gpu ID %d : Total Memory: %12.6f\n",gpu_id,double(GetAvailableMemory())*1e-9);
	


	//Some constants we might need
	MaxDegen = pmax_degen;
	SymNrepres = sym_nrepres;
	Jmax = jmax;

	//Copy constants to the GPU
	//cudaMemcpyToSymbol (g_maxdegen,&MaxDegen, sizeof(int) );		
	//cudaMemcpyToSymbol (c_jmax,&Jmax , sizeof(int) );		
	// cudaMemcpyToSymbol ( g_symnrepres,&SymNrepres, sizeof(int) );	
	CheckCudaError("Copy Constants");


	//Allocate streams
	//Half_linestrength
	for(int i = 0; i < MAX_STREAMS; i++){
		half_ls_stream.push_back(NULL);
		cudaStreamCreate(&half_ls_stream.back());
	}
	

	//DotProduct
	for(int i =0; i < Nprocs; i++){
		dot_product_omp_stream.push_back(NULL);
		cudaStreamCreate(&dot_product_omp_stream.back());
	}

	Log("Transfering ThreeJ symbols\n");
	//Transfer threeJ symbols
	
	



	Log("..Done!\n");


}

void GpuManager::TransferBasisSet(const BasisSet* basisSet){
	cudaSetDevice(gpu_id);
	GpuBasisSet tmp_bset;
	
	tmp_bset.J = basisSet->GetJval();
	
	Log("Transferring Basis Set J=%d\n",basisSet->GetJval());

	tmp_bset.dimensions = basisSet->GetDimensions();

	//K Tau
	AllocateGpuMemory(&tmp_bset.KTau,sizeof(int)*size_t(tmp_bset.dimensions));
	TransferToGpu(tmp_bset.KTau, basisSet->GetKTau());

	//KBlocks
	AllocateGpuMemory(&tmp_bset.Kblock_size,sizeof(int)*size_t(tmp_bset.J+1));
	TransferToGpu(tmp_bset.Kblock_size, basisSet->GetKBlock());	
	//Vibrational index
	AllocateGpuMemory(&tmp_bset.vib_index,sizeof(int)*size_t(tmp_bset.dimensions));
	TransferToGpu(tmp_bset.vib_index, basisSet->GetVibIndex());	

	//We need the starting K information here
	tmp_bset.KStart.assign(basisSet->GetKStart(),basisSet->GetKStart()+tmp_bset.J +1);

	//Now lets put everything nicely here
	basisSets.push_back(tmp_bset);
	
	


}

void AllocateVectors(int nJ,int nsizemax,int dimenmax){
	cudaSetDevice(gpu_id);
	Log("Allocating Vectors in GPU %d...",gpu_id);
	NsizeMax = nsizemax;
	DimenMax = dimenMax;	
	
	//Allocate Host memory

	//Allocate pinned memory	
	AllocatePinnedMemory(&host_vectorI,sizeof(double)*size_t(NsizeMax));
	//allocate for each Nproc
	for(int i = 0; i <Nprocs; i++){
		host_vectorF.push_back(NULL);
		AllocatePinnedMemory(&host_vectorF.back(),sizeof(double)*size_t(NsizeMax));
	}
	
	//allocate host linstrenght
	for(int i = 0; i <Nprocs; i++){
		host_linestrength.push_back(NULL);
		AllocatePinnedMemory(&host_linestrength.back(),sizeof(double)*size_t(MaxDegen*MaxDegen));
	}



	//allocate GPU memory

	/////////////////////////Half-linestrength vectors//////////////////////
	AllocateGpuMemory(&vectorI,sizeof(double)*size_t(NsizeMax));

	//Allocate for each J and each degeneracy
	for(int i =0; i < nJ; i++){
		prim_half_ls_vectors.push_back(std::vector<double*>());
		half_ls_vectors.push_back(std::vector<double*>());
		for(int deg = 0; deg < MaxDeg; deg++){
			prim_half_ls_vectors.back().push_back(NULL);
			half_ls_vectors.back().push_back(NULL);
			AllocateGpuMemory(&prim_half_ls_vectors.back().back(),sizeof(double)*size_t(DimenMax));
			AllocateGpuMemory(&half_ls_vectors.back().back(),sizeof(double)*size_t(DimenMax));
		}
	}

	

	/////////////////final vectors and linestrength
	for(int i = 0; i <Nprocs; i++){
		vectorF.push_back(NULL);
		AllocateGpuMemory(&vectorF.back(),sizeof(double)*size_t(NsizeMax));
		prim_vectorF.push_back(NULL);
		AllocateGpuMemory(&prim_vectorF.back(),sizeof(double)*size_t(DimenMax));
		linestrength.push_back(NULL);
		AllocateGpuMemory(&linestrength.back(),sizeof(double)*size_t(MaxDegen*MaxDegen));
	}
	Log("Done!!!");
}


void GpuManager::TransferInflation(int* icontr_, int* ijterm,int dimen,std::vector<int> Ntot,std::vector<int> sym_degen){
	cudaSetDevice(gpu_id);
	GpuInflation tmp_inflate;
	
	


}
void GpuManager::TransferDipole(Dipole* dipole_){
	cudaSetDevice(gpu_id);

	dipole_me = dipole_;

	//Initialize it
	//dipole


}

	

void GpuManager::ExecuteHalfLs(double* half_ls,double * vector, int N,int indI,int indF,int idegI,int jF,int igammaI){
	cudaSetDevice(gpu_id);

}
void GpuManager::ExecuteDotProduct(double* int indI,int indF,int idegI,int jF,int igammaI){
	cudaSetDevice(gpu_id);
}
