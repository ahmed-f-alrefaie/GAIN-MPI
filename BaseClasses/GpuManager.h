
#include "../common/BaseProcess.h"
#include "../common/BaseManager.h"
#include "../common/defines.h"
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
	std::vector<double*> repres;
	
};


class GpuManager : public BaseProcess ,public BaseManager {
protected:
	int gpu_id;

	int nprocs

	Dipole* dipole_me;

	int block_in_mem;

	//Basis set information
	std::vector<GpuBasisSet> basisSets;
	std::vector<GpuInflation> inflationData;	
	
	//ThreejSymbos
	double* threejsymbols;

	//////////////////Host memory///////////////////
	double* host_vectorI;
	std:vector<double*> host_vectorF;
	std::vector<double*> host_linestrength;


	/////////////////Gpu Memory//////////////////////////
	double* vectorI;
	std::vector< std:vector<double*> > prim_half_ls_vectors;
	std::vector< std:vector<double*> > half_ls_vector;
	//Our primitive expanded vectors for each omp_thread
	std:vector<double*> prim_vectorF;
	std::vector<double*> linestrength;



	
	int NsizeMax;
	int DimenMax;
	int MaxDegen;
	int SymNrepres;
	int Jmax;

	int Nprocs;	

	std::vector<cudaStream_t>  half_ls_stream;
	std::vector<cudaStream_t> dot_product_omp_stream;
	void AllocateGpuMemory(void** mem, size_t size);
	void AllocatePinnedMemory(void** mem,size_t size);
	void TransferToGpu(void* dst,void* src,size_t size);
	void TransferToHost(void* dst,void* src,size_t size);

public :
	GpuManager(int pgpu_id,int nprocs);
	void InitializeAndTransferConstants(int jmax,int sym_nrepres);
	void TransferBasisSet(const BasisSet* basisSet);
	void TransferInflation(int* icontr_, int* ijterm);
	void TransferWigner();
	void TransferDipole(Dipole* dipole_);
	void AllocateVectors(int nsizemax,int dimenmax);

	static void PinVectorMemory(double* vector,int* n);
	static void UnpinVectorMemory(double* vector,int* n);
	
	double* GetInitialVector(){return host_vectorI;};
	double* GetFinalVector(int proc_id){return host_vectorF.at(proc_id);};
	

	void ExecuteHalfLs(double* half_ls,double * vector, int N,int indI,int indF,int idegI,int jF,int igammaI);
	void ExecuteDotProduct(double* vector,int* n int indI,int indF,int idegI,int jF,int igammaI);


};




