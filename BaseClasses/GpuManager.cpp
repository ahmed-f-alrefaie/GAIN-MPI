#include "GpuManager.h"

extern "C" double c_three_j(int * j1,int * j2,int * j3,int * k1,int * k2,int * k3);

GpuManager::GpuManager(int pgpu_id,int nprocs) : BaseProcess(), BaseManager(), gpu_id(pgpu_id), Nprocs(nprocs){



}

void GpuManager::AllocateGpuMemory(void** mem, size_t size){
		if(cudaSuccess != cudaMalloc(mem,size)){
			CheckCudaError("Memory Allocation");
		}
		TrackMemory(size);
}

void GpuManager::AllocatePinnedMemory(void** mem, size_t size){
		if(cudaSuccess != cudaMallocHost(mem,size)){
			CheckCudaError("Memory Allocation");
		}
		BaseManager::TrackGlobalMemory(size);
}

void GpuManager::TransferToGpu(void* dst,const void* src,size_t size){
		cudaSetDevice(gpu_id);
		if(cudaSuccess != cudaMemcpy(dst,src,size,cudaMemcpyHostToDevice)){
			CheckCudaError("Memory Transfer H -> D");
		}
}

void GpuManager::TransferToHost(void* dst,const void* src,size_t size){
		cudaSetDevice(gpu_id);
		if(cudaSuccess != cudaMemcpy(dst,src,size,cudaMemcpyDeviceToHost)){
			CheckCudaError("Memory Transfer D -> H");
		}
}


void GpuManager::InitializeAndTransferConstants(int jmax,int sym_repres,int pmax_degen){
	//set the device

	Log("Initializing GPU device...",gpu_id);

	cudaSetDevice(gpu_id);
	cudaFree(0);
	CheckCudaError("Wake Device");
	

	//Lets get some data goin:
	cudaDeviceProp devProp;
	cudaGetDeviceProperties(&devProp, gpu_id);
	
	//Lets get some gpu_info.
	InitializeMemory(size_t(double(devProp.totalGlobalMem)*0.95));
	Log("Gpu ID %d : Total Memory: %12.6f\n",gpu_id,double(GetAvailableMemory())*1e-9);
	


	//Some constants we might need
	MaxDegen = pmax_degen;
	SymNrepres = sym_repres;
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

	Log("Transfering ThreeJ symbols");
	//Transfer threeJ symbols
	double* tmp_three_J  = new double[(jmax+1)*(jmax+1)*3*3];

        for(int jI=0; jI <= jmax; jI++)
		for(int jF=std::max(jI-1,0); jF <= std::min(jI+1,jmax); jF++)
			for(int kI=0; kI <= jI; kI++)
				for(int kF=std::max(kI-1,0); kF <= std::min(kI+1,jF); kF++)
				{
					int dum_one = 1;
					int dum_diff = kF-kI;
					int dum_nK = -kF;
					//printf("%i %i %i %i\n",jI,kI,jF-jI,kF-kI);
					//printf("Array position = %i\n",jI + kI*(jmax+1) +(jF-jI + 1)*(jmax+1)*(jmax+1) +  (kF-kI + 1)*(jmax+1)*(jmax+1)*3);
					//threej(jI, kI, jF - jI, kF - kI) = three_j(jI, 1, jF, kI, kF - kI, -kF)
					double three = c_three_j(&jI, &dum_one , &jF, &kI, &dum_diff, &dum_nK); 
					//printf("(%i,%i,%i,%i) = %14.3E\n",jI,jF,kI,kF,three);
					tmp_three_J[jI + kI*(jmax+1) +(jF-jI + 1)*(jmax+1)*(jmax+1) +  (kF-kI + 1)*(jmax+1)*(jmax+1)*3] = three;
					//printf("three-j[%i,%i,%i,%i] = %12.6f\n",jI,jF,kI,kF,three_j(jI, 1, jF, kI, kF - kI, -kF));
				}	
	
	AllocateGpuMemory((void**)&threejsymbols,sizeof(double)*size_t((jmax+1)*(jmax+1)*3*3));
	TransferToGpu(threejsymbols, tmp_three_J,sizeof(double)*size_t((jmax+1)*(jmax+1)*3*3));

	//Transfer symmetry and threeJ constants

	copy_symmetry_constants(SymNrepres,MaxDegen);
	copy_jmax_constant(Jmax);


	CheckCudaError("Initialization");

	delete [] tmp_three_J;

	Log("..Done!\n");


}

void GpuManager::TransferBasisSet(BasisSet* basisSet){
	cudaSetDevice(gpu_id);
	GpuBasisSet tmp_bset;
	
	tmp_bset.J = basisSet->GetJval();
	
	Log("Transferring Basis Set J=%d\n",basisSet->GetJval());

	tmp_bset.dimensions = basisSet->GetDimensions();

	Log("Basis set is dimension size %d\n",tmp_bset.dimensions);


	//K Tau
	AllocateGpuMemory((void**)&tmp_bset.KTau,sizeof(int)*size_t(tmp_bset.dimensions));
	TransferToGpu(tmp_bset.KTau, basisSet->GetKTau(),sizeof(int)*size_t(tmp_bset.dimensions));

	//KBlocks
	AllocateGpuMemory((void**)&tmp_bset.Kblock_size,sizeof(int)*size_t(tmp_bset.J+1));
	TransferToGpu(tmp_bset.Kblock_size, basisSet->GetKBlock(),sizeof(int)*size_t(tmp_bset.J+1));	
	//Vibrational index
	AllocateGpuMemory((void**)&tmp_bset.vib_index,sizeof(int)*size_t(tmp_bset.dimensions));
	TransferToGpu(tmp_bset.vib_index, basisSet->GetVibIndex(),sizeof(int)*size_t(tmp_bset.dimensions));	

	//We need the starting K information here
	tmp_bset.KStart.assign(basisSet->GetKStart(),basisSet->GetKStart()+tmp_bset.J +1);

	//Now lets put everything nicely here
	basisSets.push_back(tmp_bset);
	
	
	CheckCudaError("Basis Set Transfer");

}

void GpuManager::AllocateVectors(int nJ,int nsizemax,int dimenmax){
	cudaSetDevice(gpu_id);
	Log("Allocating Vectors in GPU %d...",gpu_id);
	NsizeMax = nsizemax;
	DimenMax = dimenmax;	
	
	//Allocate Host memory

	//Allocate pinned memory	
	AllocatePinnedMemory((void**)&host_vectorI,sizeof(double)*size_t(NsizeMax));
	//allocate for each Nproc
	for(int i = 0; i <Nprocs; i++){
		host_vectorF.push_back(NULL);
		AllocatePinnedMemory((void**)&host_vectorF.back(),sizeof(double)*size_t(NsizeMax));
	}
	
	//allocate host linstrenght
	for(int i = 0; i <Nprocs; i++){
		host_linestrength.push_back(NULL);
		AllocatePinnedMemory((void**)&host_linestrength.back(),sizeof(double)*size_t(MaxDegen*MaxDegen));
	}



	//allocate GPU memory

	/////////////////////////Half-linestrength vectors//////////////////////
	AllocateGpuMemory((void**)&vectorI,sizeof(double)*size_t(NsizeMax));

	//Allocate for each J and each degeneracy
	for(int i =0; i < nJ; i++){
		prim_half_ls_vectors.push_back(std::vector<double*>());
		half_ls_vectors.push_back(std::vector<double*>());
		for(int deg = 0; deg < MaxDegen; deg++){
			prim_half_ls_vectors.back().push_back(NULL);
			half_ls_vectors.back().push_back(NULL);
			AllocateGpuMemory((void**)&prim_half_ls_vectors.back().back(),sizeof(double)*size_t(DimenMax));
			AllocateGpuMemory((void**)&half_ls_vectors.back().back(),sizeof(double)*size_t(DimenMax));
		}
	}

	

	/////////////////final vectors and linestrength
	for(int i = 0; i <Nprocs; i++){
		vectorF.push_back(NULL);
		AllocateGpuMemory((void**)&vectorF.back(),sizeof(double)*size_t(NsizeMax));
		prim_vectorF.push_back(NULL);
		AllocateGpuMemory((void**)&prim_vectorF.back(),sizeof(double)*size_t(DimenMax));
		linestrength.push_back(NULL);
		AllocateGpuMemory((void**)&linestrength.back(),sizeof(double)*size_t(MaxDegen*MaxDegen));
	}
	Log("Done!!!");
}


void GpuManager::TransferInflation(int* icontr_, int* ijterm,int dimen,int maxsymcoeffs,int matsize,std::vector<double*> repres,std::vector<int*> N,std::vector<int> Ntot,std::vector<int> sym_degen){
	cudaSetDevice(gpu_id);
	GpuInflation tmp_inflate;

	Log("Transferring Inflation\n");
	tmp_inflate.MaxSymCoeffs = maxsymcoeffs;
	tmp_inflate.Ntotal = Ntot;
	
	//Transfer ijtems and icontr
	AllocateGpuMemory((void**)&tmp_inflate.ijTerms,sizeof(int)*tmp_inflate.MaxSymCoeffs*SymNrepres);
	TransferToGpu(tmp_inflate.ijTerms,ijterm,sizeof(int)*tmp_inflate.MaxSymCoeffs*SymNrepres);

	AllocateGpuMemory((void**)&tmp_inflate.contr,sizeof(int)*size_t(dimen*2));
	TransferToGpu(tmp_inflate.contr,icontr_,sizeof(int)*size_t(dimen*2));


	for(int i = 0; i < SymNrepres; i++){
		tmp_inflate.sDeg.push_back(sym_degen[i]);
		
		tmp_inflate.N.push_back(NULL);
		AllocateGpuMemory((void**)&tmp_inflate.N.back(),sizeof(int)*tmp_inflate.MaxSymCoeffs);	
		TransferToGpu(tmp_inflate.N.back(),N[i],sizeof(int)*tmp_inflate.MaxSymCoeffs);

		

		tmp_inflate.repres.push_back(NULL);
		AllocateGpuMemory((void**)&tmp_inflate.repres.back(),sizeof(double)*tmp_inflate.sDeg[i]*tmp_inflate.Ntotal[i]*matsize);	
		TransferToGpu(tmp_inflate.repres.back(),repres[i],sizeof(double)*tmp_inflate.sDeg[i]*tmp_inflate.Ntotal[i]*matsize);
	}
	
	Log("done!\n");
	
	CheckCudaError("Transfer Inflation");


}
void GpuManager::TransferDipole(Dipole* dipole_,int block){
	cudaSetDevice(gpu_id);

	Log("Transfering dipole block %d to GPU %d\n",block,gpu_id);

	dipole_me = dipole_;
	dipole_block = block;
	//This will be initialized in the MultiGpuManager Later
	//Get the availablememory
	
	//dipole_me->Initialize(GetAvailableMemory());

	//Copy dipole information anyway
	copy_dipole_constant(dipole_me->GetMaxContracts());


	if(dipole_me->IsBlocked()){
		Log("Blocking not yet implemented here!\n");
		exit(0);
	}

	//Alloc to the GPU the biggest block size
	
	AllocateGpuMemory((void**)&gpu_dipole,dipole_me->GetBiggestBlockSizeBytes());		
	//Transfer tipole
	TransferToGpu(gpu_dipole,dipole_me->GetDipolePiece(dipole_block),dipole_me->GetDipoleSizeBytes(dipole_block));


	CheckCudaError("Transfer Dipole");
	Log("Done\n");
	
}

	

void GpuManager::ExecuteHalfLs(double* half_ls,double * vector, int N,int indI,int indF,int idegI,int jF,int igammaI){
	cudaSetDevice(gpu_id);

}
void GpuManager::ExecuteDotProduct(double* vector,int N, int indI,int indF,int idegI,int jF,int igammaI,int proc,double* ls){
	cudaSetDevice(gpu_id);
}
