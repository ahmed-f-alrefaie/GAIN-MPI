
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
	int dimensions; 
	std::vector<int> KStart;
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

	//Basis set information
	std::vector<GpuBasisSet> basisSets;
	std::vector<GpuInflation> inflationData;	
	


	//////////////////Host memory///////////////////
	double* host_vectorI;
	std::vector<double*> host_vectorF;
	std::vector<double*> host_linestrength;


	/////////////////Gpu Memory//////////////////////////
	double* vectorI;
	std::vector< std::vector<double*> > prim_half_ls_vectors;
	std::vector< std::vector<double*> > half_ls_vectors;
	//Our primitive expanded vectors for each omp_thread
	std::vector<double*> vectorF;
	std::vector<double*> prim_vectorF;
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
	std::vector<cudaStream_t>  half_ls_stream;
	std::vector<cudaStream_t> dot_product_omp_stream;

	//whetherthe dipole has transfered
	cudaEvent_t transfer_dipole_event;
	//Whether the half_linestrength  has completed for a particular block
	cudaEvent_t half_ls_piece_event;
	

	void AllocateGpuMemory(void** mem, size_t size);
	void AllocatePinnedMemory(void** mem,size_t size);
	void TransferToGpu(void* dst,const void* src,size_t size);
	void TransferToHost(void* dst,const void* src,size_t size);


	//Mic variable
	
	int hls_stream_id;
	int GetHStreamId(){if(hls_stream_id>=MAX_STREAMS) hls_stream_id=0; return hls_stream_id++;};


public :
	GpuManager(int pgpu_id,int nprocs);
	void InitializeAndTransferConstants(int jmax,int sym_repres,int pmax_degen);
	void TransferBasisSet(BasisSet* basisSet);
	void TransferInflation(int* icontr_, int* ijterm,int dimen,int maxsymcoeffs,int matsize,std::vector<double*> repres,std::vector<int*> N,std::vector<int> Ntot,std::vector<int> sym_degen);
	void TransferWigner();
	void TransferDipole(Dipole* dipole_,int block);
	void AllocateVectors(int nJ,int nsizemax,int dimenmax);

	static void PinVectorMemory(double* vector,int* n);
	static void UnpinVectorMemory(double* vector,int* n);
	
	double* GetInitialVector(){return host_vectorI;};
	double* GetFinalVector(int proc_id){return host_vectorF.at(proc_id);};
	

	void ExecuteHalfLs(double* half_ls,double * vector, int N,int indI,int indF,int idegI,int jF,int igammaI);
	void ExecuteDotProduct(double* vector,int N, int indI,int indF,int idegI,int jF,int igammaI,int proc,double* ls);

	void WaitForDevice(){cudaSetDevice(gpu_id); cudaDeviceSynchronize();};

};




