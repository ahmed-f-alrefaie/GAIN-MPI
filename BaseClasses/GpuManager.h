
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
	int K;
	int* KTau;
	int* vib_index;
	int* Kblock_size;
	int dimensions; 
};

struct GpuInflation{
	int* ijTerms;
	int * contr;
	std::vector<int*> N;
	int sDeg;
	
};


class GpuManager : public BaseProcess ,public BaseManager {
protected:
	int gpu_id;

	int nprocs

	Dipole* dipole_me;
	int dipole_block	

	

	//Basis set information
	std::vector<GpuBasisSet> basisSets;
	std::vector<GpuInflation> inflationData;	
	
	double* host_vectorI;
	//Our contracted vectors for each omp_thread
	std:vector<double*> host_vectorF;
	//Our primitive expanded vectors for each omp_thread
	std:vector<double*> prim_vectorsF;

	double* vectorI;
	
	std:vector<double*> prim_half_ls_vectors;
	std:vector<double*> half_ls_vector;
	
	int NsizeMax;
	int DimenMax;
	int MaxDegen;
	int SymNrepres;
	int Jmax;

	int Nprocs;	

	std::vector<cudaStream_t>  half_ls_stream;
	std::vector<cudaStream_t> dot_product_omp_stream;

public :
	GpuManager(int pgpu_id,int nprocs);
	void InitializeAndTransferConstants(int jmax,int sym_nrepres);
	void TransferBasisSet(const BasisSet* basisSet);
	void TransferInflation(int* icontr_, int* ijterm);
	void TransferDipole(Dipole* dipole_);
	void AllocateVectors(int nsizemax,int dimenmax);

	static void PinVectorMemory(double* vector,int* n);
	static void UnpinVectorMemory(double* vector,int* n);
	
	double* GetInitialVector(){return host_

	void ExecuteHalfLs(double* half_ls,double * vector, int N,int indI,int indF,int idegI,int jF,int igammaI);
	void ExecuteDotProduct(double* vector,int* n int indI,int indF,int idegI,int jF,int igammaI);


};




