#include "GpuManager.h"
#include "../common/BaseProcess.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include "../GPU/dipole_GPU.h"
#include "Input.h"
#include "States.h"
#pragma once


class MultiGpuManager : public BaseProcess {

private:
	int total_gpus_assigned;
	std::vector<GpuManager*> m_gpus;
	Dipole* m_dipole;
	States* m_states;
	std::vector<int> m_jvals;
	int device_num;
	std::vector<int> m_dipole_dist;
	
	std::vector<int> proc_distrub;
	std::vector<int> total_procs;
	int gpu;
	//Used to get the linestrength result
	std::vector< std::vector< double* > > half_linestrength;
	//Used to store the temporary result when it is completed with blocking
	std::vector< std::vector< double* > > tmp_half_linestrength;

	int nJ;
	double* vectorI;

	int n_procs;

	int NsizeMax;
	int DimenMax;
	int DegenMax;
	//Methods
	int GetFreeDevice();

	
public:
	MultiGpuManager(std::vector<int> jvals,States* states,int nprocs,bool rotsym=false);
	
	//Transfer functions
	void InitializeAndTransferConstants(int jmax,int sym_repres,int pmax_degen);
	void TransferBasisSet(BasisSet* basisSet);
	void TransferInflation(int* icontr_, int* ijterm,int dimen,int maxsymcoeffs,int matsize,std::vector<double*> repres,std::vector<int*> N,std::vector<int> Ntot,std::vector<int> sym_degen);
	void TransferDipole(Dipole* dipole_);
	void TransferWigner(std::vector<Wigner> p_wigner);

	void AllocateVectors(int nJ,int nsizemax,int dimenmax);
	
	void UpdateHalfLinestrength(double* half_ls,int jInd,int ideg);

	double* GetHalfLineStrength(int indF,int idegI);
	
	double* GetInitialVector();
	double* GetFinalVector(int proc_id);
	double* GetLinestrength(int proc_id);
	static void PinVectorMemory(double* vector,int n);
	static void UnpinVectorMemory(double* vector,int n);
	//update eigenvectors
	void UpdateEigenVector();
	void UpdateEigenVector(int proc_id);
	
	//Do work
	void ExecuteHalfLs(int iLevelI,int indI,int indF,int idegI,int igammaI);
	void ExecuteDotProduct(int indF,int idegI,int idegF,int igammaF,int proc);
	void WaitForLineStrengthResult(int proc_id);

};


