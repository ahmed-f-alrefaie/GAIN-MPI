
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
	int* N;
	int sDeg;
	int 
};


class GpuManager : public BaseProcess ,public BaseManager {
protected:
	int gpu_id;

	int nprocs

	Dipole* dipole_me;
	int dipole_block
	
	int max_degeneracy;	

	//Basis set information
	std::vector<GpuBasisSet> basisSets;
	std::vector<GpuInflation> inflationData;	
	
	//Our contracted vectors for each omp_thread
	std:vector<double*> vectorF;

	//Our primitive expanded vectors for each omp_thread
	std:vector<double*> prim_vectorsF;

	double* vectorI
	std:vector<double*> prim_half_ls_vectors;
	std:vector<double*> half_ls_vector;
	
	int NsizeMax;
	int DimenMax;
	int max_degen;
	int sym_nrepres;

	int Nprocs;	

	cudaStream_t half_ls_stream[MAX_STREAMS];
	std::vector<cudaStream_t> dot_product_omp_stream;

public :
	GpuManager(int pgpu_id,int pmax_degen,int nprocs);
	void InitializeAndTransferConstants(int jmax,int sym_nrepres);
	void TransferBasisSet(const BasisSet* basisSet);
	void TransferInflation(int* icontr_, int* ijterm);
	void TransferDipole(Dipole* dipole_);
	void AllocateVectors(int nsizemax,int dimenmax);

	

	void ExecuteHalfLs(double* half_ls,double * vector, int N,int indI,int indF,int idegI,int jF,int igammaI);
	void ExecuteDotProduct(double* int indI,int indF,int idegI,int jF,int igammaI);


};




