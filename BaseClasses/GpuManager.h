#include "../common/Wigner.h"
#include "../common/BaseProcess.h"
#include "../common/BaseManager.h"
#include "../common/defines.h"
#include "../GPU/dipole_GPU.h"
#include <vector>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <omp.h>
#include "Dipole.h"
#include "BasisSet.h"

#pragma once

struct GpuBasisSet{
	int J; 
	int* KTau;
	int* vib_index;
	int* Kblock_size;
        int* normal_K;
	//int* normal_Ktau;
	int dimensions; 
	std::vector<int> KStart;
	std::vector<int> host_KBlock;
};

struct GpuInflation{
	int* ijTerms;
	int * contr;
	std::vector<int*> N;
	std::vector<int> sDeg;
	std::vector<int> Ntotal;
	std::vector<double*> repres;
	int MaxSymCoeffs;
	
};


class GpuManager : public BaseProcess ,public BaseManager {
private:
	int gpu_id;

	int nprocs;

	Dipole* dipole_me;

	int block_in_mem;
	int dipole_block;
	bool doing_half_ls;
	//Basis set information
	std::vector<GpuBasisSet> basisSets;
	std::vector<GpuInflation> inflationData;	
	


	//////////////////Host memory///////////////////
	double* host_vectorI;
	std::vector< std::vector<double*> > host_half_ls_vectors;
	std::vector<double*> host_vectorF;
	std::vector<double*> host_linestrength;


	/////////////////Gpu Memory//////////////////////////
	double* vectorI;
	std::vector<int*> old_roots;

	std::vector<double*> eigenvect;
	
	std::vector< std::vector<Wigner> > gpu_Wigner;


	std::vector< std::vector<double*> > prim_half_ls_vectors;
	std::vector< std::vector<double*> > unsorted_half_ls_vectors;
	std::vector< std::vector<double*> > half_ls_vectors;
	//Our primitive expanded vectors for each omp_thread
	std::vector<double*> vectorF;
	std::vector<double*> prim_vectorF;
	std::vector<double*> unsorted_vectorF;
	std::vector<double*> linestrength;

	//ThreejSymbos
	double* threejsymbols;
	double* gpu_dipole;

	
	int NsizeMax;
	int DimenMax;
	int MaxDegen;
	int SymNrepres;
	int Jmax;

	int Nprocs;	
	//Streams and events

	std::vector<cudaStream_t> dot_product_omp_stream;

	std::vector<std::vector<int> > stream_id;

	std::vector<std::vector<std::vector<cudaStream_t> > > half_ls_stream;

	std::vector<std::vector<cudaEvent_t> > correlated_event;

	std::vector<std::vector<std::vector<cudaEvent_t> > > completed_half_ls_event;
	std::vector<std::vector<cudaStream_t>  > transfer_half_ls_stream;
	std::vector<std::vector<cudaStream_t>  >correlate_half_ls_stream;

	//whetherthe dipole has transfered
	cudaEvent_t transfer_dipole_event;
	cudaStream_t transfer_dipole_stream;
	//Whether the half_linestrength  has completed for a particular block
	cudaEvent_t half_ls_piece_event;
	

	void AllocateGpuMemory(void** mem, size_t size);
	void AllocatePinnedMemory(void** mem,size_t size);
	void TransferToGpu(void* dst,const void* src,size_t size);
	void TransferToHost(void* dst,const void* src,size_t size);

	


	//Mic variable
	
	int hls_stream_id;
	int GetHStreamId(int indF,int degI){if(stream_id[indF][degI]>=MAX_STREAMS) stream_id[indF][degI]=0; return stream_id[indF][degI]++;};
	void CheckCudaError(const char* tag);

	//CublasRelated
	std::vector<cublasHandle_t> handle;
	cublasStatus_t stat;
	bool rotsym_do;
	void ExecuteKBlockHalfLs(int indI,int indF,int idegI,int igammaI);
	void ExecuteRotSymHalfLs(int indI,int indF,int idegI,int igammaI);
	
public :
	GpuManager(int pgpu_id,int nprocs,bool rotsym=false);


	//Initialization and allocation methods
	void InitializeAndTransferConstants(int jmax,int sym_repres,int pmax_degen);
	void TransferBasisSet(BasisSet* basisSet);
	void TransferInflation(int* icontr_, int* ijterm,int dimen,int maxsymcoeffs,int matsize,std::vector<double*> repres,std::vector<int*> N,std::vector<int> Ntot,std::vector<int> sym_degen);
	void TransferWigner(std::vector<Wigner> p_wigner);
	void TransferDipole(Dipole* dipole_,int block);
	void SwitchDipoleBlock(int block);
	void AllocateVectors(int nJ,int nsizemax,int dimenmax);
	



	//Get Methods
	double* GetInitialVector(){return host_vectorI;};
	double* GetFinalVector(int proc_id){return host_vectorF.at(proc_id);};
	double* GetLinestrength(int proc_id){return host_linestrength.at(proc_id);};
	void GetHalfLineStrengthResult(double* half_ls,int indF,int idegI);
	int GetCurrentBlock(){return dipole_block;};
	bool WillDoHalfLs(){return doing_half_ls;}

	//Update GPU methods
	void UpdateHalfLinestrength(double* half_ls,int jInd,int ideg);
	void UpdateEigenVector();
	void UpdateEigenVector(int proc_id);
	
	
	//Work methods
	void TransformHalfLsVector(int indI,int indF,int idegI,int igammaI);
	void ExecuteHalfLs(int indI,int indF,int idegI,int igammaI);
	void ExecuteDotProduct(int indF,int idegI,int idegF,int igammaF,int proc);

	//Wait methods
	void WaitForLineStrengthResult(int proc_id);
	void WaitForDevice(){cudaSetDevice(gpu_id); cudaDeviceSynchronize();};

};




