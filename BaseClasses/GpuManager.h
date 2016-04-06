
#include "../common/BaseProcess.h"
#include "../common/BaseManager.h"
#include "../common/defines.h"
#include <vector>
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <omp.h>
#include "../TROVE/Dipole.h"
#include "../GPU/dipole_kernals.cuh"
#pragma once


class GpuManager : public BaseProcess ,public BaseManager {
protected:
	int gpu_id;

	Dipole* dipole_me;
	int dipole_blcok

	

	//Basis set information
	std::vector<int> J; 
	std::vector<int*> K;
	std::vector<int*> KTau;
	std::vector<int*> vib_index;
	std::vector<int*> Kblock_size;
	std::vector<int> dimensions; 		
	
	//Our contracted vectors for each omp_thread
	std:vector<double*> vectors;

	//Our primitive expanded vectors for each omp_thread
	std:vector<double*> prim_vectors;

	std:vector<double*> half_ls_vector;
	



	std::vector<cudaStream_t> half_ls_stream;
	int max_degeneracy;
	std:vector<cudaStream_t> dot_product_omp_stream;

public :
	GpuManager(int gpu_id,int max_degen);
	
	void TransferBasisSet(int J_,int* K_, int KTau_,int* vib_index_, int* k_block_);
	void TransferDipole(Dipole* dipole_,int block_num);
	//void SetMaxVectorSize(int num_elems);

	

	virtual void ExecuteHalfLs(double* vector,int N,int indI,int indF,int idegI,int jF,int igammaI)=0;
	virtual void ExecuteDotProduct(int indI,int indF,int idegI,int jF,int igammaI)=0;


};




