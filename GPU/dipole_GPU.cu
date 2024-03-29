#include "dipole_GPU.h"
#include "cublas_v2.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cmath>
#include <cstdio>
#include <vector>

//Math constants
//Kepler definitions

#define CORRELATE_BLOCK_SIZE 512
#define DIPOLE_BLOCK_SIZE 128



#define WARP_SIZE 32
#define PI 3.14159265359
#define AVOGNO 6.0221415e+23
#define PLANCK 6.62606896e-27
#define VELLGT 2.99792458e+10
#define BOLTZ 1.380658e-16
#define ACOEF 3.13618923e-7
#define SQRT2 0.70710678118

//#define DEBUG




#ifdef CUBLAS_API_H_
// cuBLAS API errors
static const char *_cudaGetErrorEnum(cublasStatus_t error)
{
    switch (error)
    {
        case CUBLAS_STATUS_SUCCESS:
            return "CUBLAS_STATUS_SUCCESS";

        case CUBLAS_STATUS_NOT_INITIALIZED:
            return "CUBLAS_STATUS_NOT_INITIALIZED";

        case CUBLAS_STATUS_ALLOC_FAILED:
            return "CUBLAS_STATUS_ALLOC_FAILED";

        case CUBLAS_STATUS_INVALID_VALUE:
            return "CUBLAS_STATUS_INVALID_VALUE";

        case CUBLAS_STATUS_ARCH_MISMATCH:
            return "CUBLAS_STATUS_ARCH_MISMATCH";

        case CUBLAS_STATUS_MAPPING_ERROR:
            return "CUBLAS_STATUS_MAPPING_ERROR";

        case CUBLAS_STATUS_EXECUTION_FAILED:
            return "CUBLAS_STATUS_EXECUTION_FAILED";

        case CUBLAS_STATUS_INTERNAL_ERROR:
            return "CUBLAS_STATUS_INTERNAL_ERROR";
    }

    return "<unknown>";
}
#endif




__constant__ int dip_stride_1;
__constant__ int dip_stride_2;
__constant__ int c_jmax;
__constant__ int g_symnrepres;
__constant__ int g_maxdegen;
int host_jmax;


__global__ void device_unsort_vector(const int dimenI,const int* N_, const double* c_vec_,double* vec_);

__global__ void device_expand_vectors(const int dimenI,const int igammaI,const int maxcoeff,const int idegI,const int sdeg,const int Ntot,const int* ijterms_, const int* icontr_,const int* N_,const double* repres_, const double* vecI_,double* vec_);
__global__ void device_sort_vector(const int dimenI,const int* N_, const double* c_vec_,double* vec_);

__global__ void device_compute_1st_half_ls_flipped_dipole_shared_nontrove(
const int dimenF,const int dimenI,const int jI,const int jF,
const int kF,const int* tauI_,const int* tauF_, const int* icorrI_,const int* icorrF_,const int startF_idx,int startI_idx,const int* kblock_size_,
const double* __restrict__ dipole_me,
const double* vector,const double*  threej,double*  half_ls);

__global__ void device_compute_1st_half_ls_flipped_dipole_shared_rotsym_nontrove(
const int dimenF,const int dimenI,const int nlevelI,const int nlevelF,const int* kI_,
const int* kF_,const int* ktauI_,const int* ktauF_, const int* icorrI_,const int* icorrF_,
const double* __restrict__ dipole_me,const double* __restrict__ wigner_,
const double* vector,double*  half_ls);

__global__ void device_compute_1st_half_ls_flipped_dipole_shared_block_nontrove(
const int dimenF,const int dimenI,const int jI,const int jF,
const int kF,const int* tauI_,const int* tauF_, const int* icorrI_,const int* icorrF_,const int startF_idx,int startI_idx,const int startFblock,const int endFblock,const int ncontrF,const int* kblock_size_,
const double* __restrict__ dipole_me,
const double* vector,const double* __restrict__ threej,double*  half_ls);


__global__ void device_compute_1st_half_ls_flipped_dipole_shared_rotsym_block_nontrove(
const int dimenF,const int dimenI,const int jI,const int jF,
const int kF,const int* tauI_,const int* tauF_, const int* icorrI_,const int* icorrF_,const double* eigenvectI_,const double* eigenvectF_,const int startF_idx,int startI_idx,const int startFblock,const int endFblock,const int ncontrF,const int* kblock_size_,
const double* __restrict__ dipole_me,
const double* vector,const double*  threej,double*  half_ls);

__global__ void device_compute_1st_half_ls_flipped_dipole_shared_rotsym_nontrove(
const int dimenF,const int dimenI,const int nlevelI,const int nlevelF,const int sym_max_degen,const int* kI_,
const int* kF_,const int* ktauI_,const int* ktauF_, const int* icorrI_,const int* icorrF_,const int startFblock,const int endFblock,const int ncontrF,
const double* __restrict__ dipole_me,const double* __restrict__ wigner_,
const double* vector,double*  half_ls);

__global__ void device_compute_1st_half_ls_flipped_dipole_shared_blocks_fast(const int dimenF,const int dimenI,const int j0dimen,const int jI,const int jF,
const int kF,const int* tauI_,const int* tauF_, const int* icorrI_,const int* icorrF_,const int startFblock,const int endFblock,const int ncontrF,
const double* __restrict__ dipole_me,
const double* vector,const double* __restrict__  threej,double*  half_ls);


void copy_symmetry_constants(int sym_nrepres,int maxdeg){
	
	cudaMemcpyToSymbol (g_symnrepres,&sym_nrepres , sizeof(int) );
	cudaMemcpyToSymbol (g_maxdegen,&maxdeg , sizeof(int) );
}

void copy_jmax_constant(int jmax){
	cudaMemcpyToSymbol (c_jmax,&jmax , sizeof(int) );

}

void copy_dipole_constant(int dip_stride){
	
	cudaMemcpyToSymbol (dip_stride_1,&dip_stride , sizeof(int) );

}


void transform_vector_primitive(const int dimenI,const int igammaI,const int maxcoeff,const int idegI,const int sdeg,const int Ntot,const int* ijterms_, const int* icontr_,const int* N_,const double* repres_, const double* vecI_,double* vec_,cudaStream_t stream){

	int gridSize = (int)ceil((float)dimenI/CORRELATE_BLOCK_SIZE);

	device_expand_vectors<<<gridSize,CORRELATE_BLOCK_SIZE,0,stream>>>(
		dimenI,
		igammaI,
		maxcoeff,
		idegI,
		sdeg,
		Ntot,
		ijterms_,
		icontr_,
		N_,
		repres_,
		vecI_,
		vec_);	


}

void sort_vector_correlated(const int dimenI,const int* N_, const double* c_vec_,double* vec_,cudaStream_t stream){

	int gridSize = (int)ceil((float)dimenI/CORRELATE_BLOCK_SIZE);

	device_sort_vector<<<gridSize,CORRELATE_BLOCK_SIZE,0,stream>>>(
		dimenI,
		N_,
		c_vec_,
		vec_);	


}

void unsort_vector_correlated(const int dimenI,const int* N_, const double* c_vec_,double* vec_,cudaStream_t stream){

	int gridSize = (int)ceil((float)dimenI/CORRELATE_BLOCK_SIZE);

	device_unsort_vector<<<gridSize,CORRELATE_BLOCK_SIZE,0,stream>>>(
		dimenI,
		N_,
		c_vec_,
		vec_);	


}


void compute_gpu_half_linestrength_(const int dimenF,const int dimenI,const int jI,const int jF,
const int kF,const int* tauI_,const int* tauF_, const int* icorrI_,const int* icorrF_,const int startF_idx,int startI_idx,const int startFblock,const int endFblock,const int ncontrF,const int kFBlocksize,const int* kblock_size_,
const double* dipole_me,
const double* vector,const double*  threej,double*  half_ls,cudaStream_t stream){


	int half_grid_size = (int)ceil((float)(dimenF)/(float)DIPOLE_BLOCK_SIZE);
				//printf("half_grid params: b:%i t:%i N:%i\n",half_grid_size,DIPOLE_BLOCK_SIZE,DIPOLE_BLOCK_SIZE*half_grid_size);
	//device_compute_1st_half_ls_flipped_dipole_shared_nontrove<<<half_grid_size,DIPOLE_BLOCK_SIZE,0,half_ls_stream[stream_id]>>>(h_k_blocks[indF][k],dimen[indI],jI,
	//										jF,k,tau[indI],tau[indF],icorr[indI],icorr[indF],
	//										k_start[indF][k],k_start[indI][max(k-1,0)],k_block_size[indI],
	//										cuda_dipole_me,
	//										corr_vector,g_threeJ,
	//										g_half_ls[indF][ideg]);
				//cudaDeviceSynchronize();
	//			CheckCudaError("K-block");
	
	
	device_compute_1st_half_ls_flipped_dipole_shared_block_nontrove<<<half_grid_size,DIPOLE_BLOCK_SIZE,0,stream>>>(
	dimenF,dimenI,jI,jF,
	kF,tauI_,tauF_,icorrI_,icorrF_,startF_idx,startI_idx,startFblock,endFblock,ncontrF,kblock_size_,
	dipole_me,
	vector,threej,half_ls);


}

void compute_gpu_half_linestrength_fast(const int dimenF,const int dimenI,const int j0dimen,const int jI,const int jF,
const int kF,const int* tauI_,const int* tauF_, const int* icorrI_,const int* icorrF_,const int startFblock,const int endFblock,const int ncontrF,
const double* dipole_me,
const double* vector,const double*   threej,double*  half_ls,cudaStream_t stream){
	
	int half_grid_size = (int)ceil((float)(dimenF)/(float)DIPOLE_BLOCK_SIZE);
				//printf("half_grid params: b:%i t:%i N:%i\n",half_grid_size,DIPOLE_BLOCK_SIZE,DIPOLE_BLOCK_SIZE*half_grid_size);
	//device_compute_1st_half_ls_flipped_dipole_shared_nontrove<<<half_grid_size,DIPOLE_BLOCK_SIZE,0,half_ls_stream[stream_id]>>>(h_k_blocks[indF][k],dimen[indI],jI,
	//										jF,k,tau[indI],tau[indF],icorr[indI],icorr[indF],
	//										k_start[indF][k],k_start[indI][max(k-1,0)],k_block_size[indI],
	//										cuda_dipole_me,
	//										corr_vector,g_threeJ,
	//										g_half_ls[indF][ideg]);
				//cudaDeviceSynchronize();
	//			CheckCudaError("K-block");
	
	
	device_compute_1st_half_ls_flipped_dipole_shared_blocks_fast<<<half_grid_size,DIPOLE_BLOCK_SIZE,0,stream>>>(
	dimenF,dimenI,j0dimen,jI,jF,
	kF,tauI_,tauF_,icorrI_,icorrF_,startFblock,endFblock,ncontrF,
	dipole_me,
	vector,threej,half_ls);


}

void compute_gpu_half_linestrength_rotsym_(
const int dimenF,const int dimenI,const int nlevelI,const int nlevelF,const int sym_max_degen,const int* kI_,
const int* kF_,const int* ktauI_,const int* ktauF_, const int* icorrI_,const int* icorrF_,const int startFblock,const int endFblock,const int ncontrF,
const double* dipole_me,const double* wigner_,
const double* vector,double*  half_ls,cudaStream_t stream){


	int half_grid_size = (int)ceil((float)(dimenF)/(float)DIPOLE_BLOCK_SIZE);
				//printf("half_grid params: b:%i t:%i N:%i\n",half_grid_size,DIPOLE_BLOCK_SIZE,DIPOLE_BLOCK_SIZE*half_grid_size);
	//device_compute_1st_half_ls_flipped_dipole_shared_nontrove<<<half_grid_size,DIPOLE_BLOCK_SIZE,0,half_ls_stream[stream_id]>>>(h_k_blocks[indF][k],dimen[indI],jI,
	//										jF,k,tau[indI],tau[indF],icorr[indI],icorr[indF],
	//										k_start[indF][k],k_start[indI][max(k-1,0)],k_block_size[indI],
	//										cuda_dipole_me,
	//										corr_vector,g_threeJ,
	//										g_half_ls[indF][ideg]);
				//cudaDeviceSynchronize();
	//			CheckCudaError("K-block");
	
	
	


	device_compute_1st_half_ls_flipped_dipole_shared_rotsym_nontrove<<<half_grid_size,DIPOLE_BLOCK_SIZE,0,stream>>>(dimenF,dimenI,nlevelI,nlevelF,sym_max_degen,kI_,kF_,ktauI_,ktauF_,icorrI_,icorrF_,startFblock,endFblock,ncontrF,dipole_me,wigner_,vector,half_ls);


}


__global__ void device_compute_1st_half_ls_flipped_dipole_shared_nontrove(
const int dimenF,const int dimenI,const int jI,const int jF,
const int kF,const int* tauI_,const int* tauF_, const int* icorrI_,const int* icorrF_,const int startF_idx,int startI_idx,const int* kblock_size_,
const double* __restrict__ dipole_me,
const double* vector,const double*  threej,double*  half_ls)
{

	//const int irootF = blockIdx.x * blockDim.x + threadIdx.x;
	//double sq2 = 1.0/sqrt(2.0);
	__shared__ double s_ls_factor[DIPOLE_BLOCK_SIZE];
	__shared__ int s_icontrI[DIPOLE_BLOCK_SIZE];
	__shared__ int s_tauI[DIPOLE_BLOCK_SIZE];
	__shared__ int s_sigmaI[DIPOLE_BLOCK_SIZE];

	int t_id = threadIdx.x;	
	//int b_size=32;

	//int b_start = (threadIdx.x/32)*32;
	//int w_id = threadIdx.x % 32;
	int g_id = blockIdx.x*blockDim.x + threadIdx.x;
	int icontrF,kI, tauI,tauF,sigmaF, sigmaI,dipole_idx,irootI,kBlockSize;
	//These are o remove if statements
	int kI_kF_eq,tauF_tauI_neq;
	double ls = 0.0,f3j=0.0,final_half_ls,sq_factor;
	int irootF=blockIdx.x*blockDim.x + threadIdx.x + startF_idx;
	bool valid;

	//if(g_id==0) printf("c_jmax:%i kF:%i startF_idx=%i endF_idx=%i blocksize=%i\n",c_jmax,kF,startF_idx,startF_idx+dimenF,dimenF);
	tauF = 0;
	final_half_ls = 0.0;
	if(irootF<startF_idx+dimenF){

		icontrF =icorrF_[irootF]-1;
		tauF  =  tauF_[irootF] & 1;
	}

	valid = irootF<(startF_idx+dimenF);

	sigmaF = (kF % 3)*tauF;
	int count_dimen;
	//Delta K = -1
	kI = max(kF-1,0);	
	for(kI=kI; kI<kF+2 && kI<(jI+1); kI++){
		//if(kblock_size_ != NULL){
		kBlockSize = kblock_size_[kI]; 
		//}
		kI*kF!=0 ? sq_factor=SQRT2 : sq_factor=1.0;

		//if(t_id==0) printf("[%i] K:%i start_idx=%i block_size=%i\n", blockIdx.x,kI,startI_idx,kBlockSize);
		kI_kF_eq = (kF==kI); 
		//kI_kF_zero = ((kI*kF) != 0);
		
		f3j  =  threej[jI + kI*(c_jmax+1) + (jF - jI + 1)*(c_jmax+1)*(c_jmax+1) + (kF - kI +1)*(c_jmax+1)*(c_jmax+1)*3];	
		//if(g_id==0) printf("threeJ[%i,%i,%i,%i]=%12.6f\n", jI,kI,jF-jI,kF-kI,f3j);
		for(int b_irootI=0; b_irootI < kBlockSize; b_irootI+=DIPOLE_BLOCK_SIZE){
			int local_idx = b_irootI+t_id;
			irootI = b_irootI+t_id+startI_idx;
			s_ls_factor[t_id]=0.0;
			s_icontrI[t_id]=-1;
				
			if(local_idx < kBlockSize){
				s_tauI[t_id] = tauI_[irootI];
				s_icontrI[t_id] = icorrI_[irootI]-1;

				s_sigmaI[t_id] = (kI % 3);

				//ls =f3j*vector[irootI]*(1.0 + (SQRT2 - 1.0)*double(kI_kF_zero)*(!kI_kF_eq));
				//printf("%i ls =%14.3e f3j=%12.6f vector=%14.3e\n",t_id,ls,f3j,vector[irootI]);

				ls = vector[irootI];		
				//Will not diverge		
				if(kF!=kI){
					ls*=sq_factor;
				}
				s_ls_factor[t_id] = ls*f3j;//*(1.0 + (SQRT2 - 1.0)*double(kI_kF_zero)*(!kI_kF_eq));//ls;
			}
			__syncthreads();
			for(int i = 0; i < DIPOLE_BLOCK_SIZE; i++){

					
				ls= s_ls_factor[i];

					
				if(ls==0.0) continue;
				//if(local_idx >= kBlockSize) continue;
				tauI=s_tauI[i];
				tauF_tauI_neq = (tauF!=tauI);
				if(valid){
				/*	if(kI==kF)
						ls*=double(tauF-tauI)*dipole_me[icontrF + s_icontrI[i]*dip_stride_1 + 2*dip_stride_2];
					else if(tauF!=tauI)
						ls*=double(kF-kI)*dipole_me[icontrF + s_icontrI[i]*dip_stride_1 + 0*dip_stride_2];
					else
						ls*=-dipole_me[icontrF + s_icontrI[i]*dip_stride_1 + 1*dip_stride_2];
				*/		 


					dipole_idx=2*(kI_kF_eq)+ (!kI_kF_eq)*(!tauF_tauI_neq);

					ls *= double((tauF-tauI)*(kI_kF_eq) + (tauF-tauI)*(kF-kI)*(!kI_kF_eq)*( tauF_tauI_neq) - (!kI_kF_eq)*(!tauF_tauI_neq));

					sigmaI = s_sigmaI[i]*tauI;
					sigmaI = 2*(~(sigmaI+kI) & 1)-1;
					final_half_ls+=ls*double(sigmaI)*dipole_me[icontrF + s_icontrI[i]*dip_stride_1 +dipole_idx*dip_stride_2];
;//f3j;//ls;// double(sigmaI)*ls*dipole_me[icontrF + s_icontrI[i + b_start]*dip_stride_1 + dipole_idx*dip_stride_2];
				}
					
			} 
			__syncthreads();

		}
		startI_idx+=kBlockSize;
		if(startI_idx>=dimenI) break;

	}
	if(valid)
	{
		final_half_ls *= double(2*(~(sigmaF) & 1)-1);
		half_ls[irootF] = final_half_ls;
	}	

	

}


__global__ void device_compute_1st_half_ls_flipped_dipole_shared_block_nontrove(
const int dimenF,const int dimenI,const int jI,const int jF,
const int kF,const int* tauI_,const int* tauF_, const int* icorrI_,const int* icorrF_,const int startF_idx,int startI_idx,const int startFblock,const int endFblock,const int ncontrF,const int* kblock_size_,
const double* __restrict__ dipole_me,
const double* vector,const double* __restrict__  threej,double*  half_ls)
{

	//shared memory to chache vectors and quanta
	__shared__ double s_ls_factor[DIPOLE_BLOCK_SIZE];
	__shared__ int s_icontrI[DIPOLE_BLOCK_SIZE];
	__shared__ int s_tauI[DIPOLE_BLOCK_SIZE];
	__shared__ double s_sigmaI[DIPOLE_BLOCK_SIZE];

	//The thread_id in the block
	int t_id = threadIdx.x;	
	//The global id in the grid
	int g_id = blockIdx.x*blockDim.x + threadIdx.x;
	int stride = ncontrF*dip_stride_1;
	//Lower and uper quantum numebrs
	int icontrF,kI, tauI,tauF,sigmaF, sigmaI,dipole_idx,irootI,kBlockSize;

	//Used to determine which dipole without branching	
	int kI_kF_eq,tauF_tauI_neq,kI_kF_neq;

	double ls = 0.0,f3j=0.0,final_half_ls,sq_factor;
	int irootF=blockIdx.x*blockDim.x + threadIdx.x + startF_idx;
	bool valid;

	tauF = 0;
	final_half_ls = 0.0;
	if(irootF<startF_idx+dimenF){

		icontrF =icorrF_[irootF];
		tauF  =  tauF_[irootF] & 1;
	}


	valid = irootF<(startF_idx+dimenF) && ((icontrF >=startFblock) && (icontrF < endFblock));
	icontrF -= startFblock;
	sigmaF = (kF % 3)*tauF;
	int count_dimen;
	//Delta K = -1
	kI = max(kF-1,0);	
	for(kI=kI; kI<kF+2 && kI<(jI+1); kI++){

		kBlockSize = kblock_size_[kI]; 

		kI*kF!=0 ? sq_factor=SQRT2 : sq_factor=1.0;


		kI_kF_eq = (kF==kI); 
		kI_kF_neq = (kF!=kI); 
		
		f3j  =  threej[jI + kI*(c_jmax+1) + (jF - jI + 1)*(c_jmax+1)*(c_jmax+1) + (kF - kI +1)*(c_jmax+1)*(c_jmax+1)*3];	

		for(int b_irootI=0; b_irootI < kBlockSize; b_irootI+=DIPOLE_BLOCK_SIZE){
			int local_idx = b_irootI+t_id;
			irootI = b_irootI+t_id+startI_idx;
			s_ls_factor[t_id]=0.0;
			s_icontrI[t_id]=-1;
				
			if(local_idx < kBlockSize){
				s_tauI[t_id] = tauI_[irootI];
				s_icontrI[t_id] = icorrI_[irootI]*ncontrF;

				sigmaI = (kI % 3)*s_tauI[t_id];

				sigmaI = 2*(~(sigmaI+kI) & 1)-1;
				s_sigmaI[t_id] = double(sigmaI);

				ls = vector[irootI];		
				//Will not diverge		
				if(kF!=kI){
					ls*=sq_factor;
				}
				s_ls_factor[t_id] = ls*f3j;//*(1.0 + (SQRT2 - 1.0)*double(kI_kF_zero)*(!kI_kF_eq));//ls;
			}
			__syncthreads();
			for(int i = 0; i < DIPOLE_BLOCK_SIZE; i++){

					
				ls= s_ls_factor[i];

					
				if(ls==0.0) continue;
				//if(local_idx >= kBlockSize) continue;
				tauI=s_tauI[i];
				tauF_tauI_neq = (tauF!=tauI);
				if(valid){
				/*	if(kI==kF)
						ls*=double(tauF-tauI)*dipole_me[icontrF + s_icontrI[i]*dip_stride_1 + 2*dip_stride_2];
					else if(tauF!=tauI)
						ls*=double(kF-kI)*dipole_me[icontrF + s_icontrI[i]*dip_stride_1 + 0*dip_stride_2];
					else
						ls*=-dipole_me[icontrF + s_icontrI[i]*dip_stride_1 + 1*dip_stride_2];
				*/		 


					dipole_idx=2*(kI_kF_eq)+ (kI_kF_neq)*(!tauF_tauI_neq);

					ls *= double((tauF-tauI)*(kI_kF_eq) + (tauF-tauI)*(kF-kI)*(kI_kF_neq)*( tauF_tauI_neq) - (kI_kF_neq)*(!tauF_tauI_neq));


					final_half_ls+=ls*s_sigmaI[i]*dipole_me[icontrF + s_icontrI[i] +dipole_idx*stride];
;//f3j;//ls;// double(sigmaI)*ls*dipole_me[icontrF + s_icontrI[i + b_start]*dip_stride_1 + dipole_idx*dip_stride_2];
				}
					
			} 
			__syncthreads();

		}
		startI_idx+=kBlockSize;
		if(startI_idx>=dimenI) break;

	}
	if(valid)
	{
		final_half_ls *= double(2*(~(sigmaF) & 1)-1);
		half_ls[irootF] = final_half_ls;
	}	

	

}

/*
__global__ void device_compute_1st_half_ls_flipped_dipole_shared_rotsym_block_nontrove(
const int dimenF,const int dimenI,const int jI,const int jF,
const int kF,const int* tauI_,const int* tauF_, const int* icorrI_,const int* icorrF_,const double* eigenvectI_,const double* eigenvectF_,const int startF_idx,int startI_idx,const int startFblock,const int endFblock,const int ncontrF,const int* kblock_size_,
const double* __restrict__ dipole_me,
const double* vector,const double*  threej,double*  half_ls)
{

	//shared memory to chache vectors and quanta
	__shared__ double s_ls_factor[DIPOLE_BLOCK_SIZE];
	__shared__ int s_icontrI[DIPOLE_BLOCK_SIZE];
	__shared__ int s_tauI[DIPOLE_BLOCK_SIZE];
	__shared__ int s_sigmaI[DIPOLE_BLOCK_SIZE];
	//The thread_id in the block
	int t_id = threadIdx.x;	
	//The global id in the grid
	int g_id = blockIdx.x*blockDim.x + threadIdx.x;

	//Lower and uper quantum numebrs
	int icontrF,kI, tauI,tauF,sigmaF, sigmaI,dipole_idx,irootI,kBlockSize;

	//Used to determine which dipole without branching	
	int kI_kF_eq,tauF_tauI_neq;

	double ls = 0.0,f3j=0.0,final_half_ls,sq_factor;
	double eigenvectF;
	int irootF=blockIdx.x*blockDim.x + threadIdx.x + startF_idx;
	bool valid;

	tauF = 0;
	final_half_ls = 0.0;
	if(irootF<startF_idx+dimenF){

		icontrF =icorrF_[irootF];
		tauF  =  tauF_[irootF] & 1;
		eigenvectF = eigenvectF_[irootF];
	}

	valid = irootF<(startF_idx+dimenF) && ((icontrF >=startFblock) && (icontrF < endFblock));

	sigmaF = (kF % 3)*tauF;
	int count_dimen;
	//Delta K = -1
	kI = max(kF-1,0);	
	for(kI=kI; kI<kF+2 && kI<(jI+1); kI++){

		kBlockSize = kblock_size_[kI]; 

		kI*kF!=0 ? sq_factor=SQRT2 : sq_factor=1.0;


		kI_kF_eq = (kF==kI); 

		
		f3j  =  threej[jI + kI*(c_jmax+1) + (jF - jI + 1)*(c_jmax+1)*(c_jmax+1) + (kF - kI +1)*(c_jmax+1)*(c_jmax+1)*3];	

		for(int b_irootI=0; b_irootI < kBlockSize; b_irootI+=DIPOLE_BLOCK_SIZE){
			int local_idx = b_irootI+t_id;
			irootI = b_irootI+t_id+startI_idx;
			s_ls_factor[t_id]=0.0;
			s_icontrI[t_id]=-1;
				
			if(local_idx < kBlockSize){
				s_tauI[t_id] = tauI_[irootI];
				s_icontrI[t_id] = icorrI_[irootI];

				s_sigmaI[t_id] = (kI % 3);
				
				//ls =f3j*vector[irootI]*(1.0 + (SQRT2 - 1.0)*double(kI_kF_zero)*(!kI_kF_eq));
				//if(g_id==0)printf("%i ls =%14.3e f3j=%12.6f vector=%14.3e\n",t_id,ls,f3j,vector[irootI]);

				ls = vector[irootI];		
				//Will not diverge		
				if(kF!=kI){
					ls*=sq_factor;
				}
				s_ls_factor[t_id] = ls*f3j*eigenvectI_[irootI];//*(1.0 + (SQRT2 - 1.0)*double(kI_kF_zero)*(!kI_kF_eq));//ls;
			}
			__syncthreads();
			for(int i = 0; i < DIPOLE_BLOCK_SIZE; i++){

					
				ls= s_ls_factor[i];

					
				if(ls==0.0) continue;
				//if(local_idx >= kBlockSize) continue;
				tauI=s_tauI[i];
				tauF_tauI_neq = (tauF!=tauI);
				if(valid){



					dipole_idx=2*(kI_kF_eq)+ (!kI_kF_eq)*(!tauF_tauI_neq);

					ls *= double((tauF-tauI)*(kI_kF_eq) + (tauF-tauI)*(kF-kI)*(!kI_kF_eq)*( tauF_tauI_neq) - (!kI_kF_eq)*(!tauF_tauI_neq));

					sigmaI = s_sigmaI[i]*tauI;
					sigmaI = 2*(~(sigmaI+kI) & 1)-1;
					final_half_ls+=ls*double(sigmaI)*dipole_me[icontrF-startFblock + s_icontrI[i]*ncontrF +dipole_idx*dip_stride_1*ncontrF]*eigenvectF;
;//f3j;//ls;// double(sigmaI)*ls*dipole_me[icontrF + s_icontrI[i + b_start]*dip_stride_1 + dipole_idx*dip_stride_2];
				}
					
			} 
			__syncthreads();

		}
		startI_idx+=kBlockSize;
		if(startI_idx>=dimenI) break;

	}
	if(valid)
	{
		final_half_ls *= double(2*(~(sigmaF) & 1)-1);
		half_ls[irootF] = final_half_ls;
	}	

	

}
*/



__global__ void device_compute_1st_half_ls_flipped_dipole_shared_rotsym_nontrove(
const int dimenF,const int dimenI,const int nlevelI,const int nlevelF,const int sym_max_degen,const int* kI_,
const int* kF_,const int* ktauI_,const int* ktauF_, const int* icorrI_,const int* icorrF_,const int startFblock,const int endFblock,const int ncontrF,
const double* __restrict__ dipole_me,const double* __restrict__ wigner_,
const double* vector,double*  half_ls)
{

	//const int irootF = blockIdx.x * blockDim.x + t0hreadIdx.x;
	//double sq2 = 1.0/sqrt(2.0);
	__shared__ double s_vector[DIPOLE_BLOCK_SIZE];
	__shared__ int s_icontrI[DIPOLE_BLOCK_SIZE];
	__shared__ int s_irlevelI[DIPOLE_BLOCK_SIZE];
	__shared__ int s_irdegI[DIPOLE_BLOCK_SIZE];

	int t_id = threadIdx.x;	
	//int b_size=32;
	double final_half_ls=0.0,dip_1,dip_2,dip_3,f_1,f_2,f_3,w_idx;
	int irootF=blockIdx.x*blockDim.x + threadIdx.x;
	int icontrI,irdegI,irlevelI,icontrF,irlevelF,irdegF,irootI;
	
	size_t dipole_idx = size_t(ncontrF)*size_t(dip_stride_1);	
	size_t wigner_idx;
	bool valid = false;
	if(irootF<dimenF){
		valid = true;
		icontrF =icorrF_[irootF];
		irlevelF  =  ktauF_[irootF];
		irdegF   = kF_[irootF];
	}

	valid &= ((icontrF >=startFblock) && (icontrF < endFblock));
	icontrF-=startFblock;
	for(int b_irootI=0; b_irootI < dimenI; b_irootI+=DIPOLE_BLOCK_SIZE){
		irootI = b_irootI+t_id;
		if(irootI < dimenI){
			s_icontrI[t_id]=icorrI_[irootI]*ncontrF;
			s_irlevelI[t_id]=ktauI_[irootI];
			s_irdegI[t_id]=kI_[irootI];
			s_vector[t_id]=vector[irootI];
		}
		__syncthreads();
		for(int i = 0; i < DIPOLE_BLOCK_SIZE; i++){

			if(i+b_irootI >= dimenI) continue;		
			if(valid){
				icontrI = s_icontrI[i];
				irlevelI= s_irlevelI[i];
				irdegI = s_irdegI[i];	
				dip_1 = dipole_me[icontrF + icontrI]; //(icontrF,icontrI,1)
				dip_2 = dipole_me[icontrF + icontrI +   dipole_idx]; //(icontrF,icontrI,2)
				dip_3 = dipole_me[icontrF + icontrI +   dipole_idx + dipole_idx]; //(icontrF,icontrI,3)
				wigner_idx = size_t(irlevelI)*3l + size_t(irlevelF)*size_t(nlevelI)*3l + size_t(irdegI)*size_t(nlevelI)*size_t(nlevelF)*3l + size_t(irdegF)*size_t(sym_max_degen)*size_t(nlevelI)*size_t(nlevelF)*3l;				
				
				//[n + ilevelI*3 + ilevelF*wigner.back().nlevelsI*3 + idegI*wigner.back().nlevelsF*wigner.back().nlevelsI*3
				//					+ idegF*sym_maxdegen*wigner.back().nlevelsF*wigner.back().nlevelsI*3]
				//f_1=0;
				//for(int n = 0; n < 3; n+){ 				
				f_1=wigner_[wigner_idx]*dip_1;
				f_2=wigner_[1+wigner_idx]*dip_2;
				f_3=wigner_[2+wigner_idx]*dip_3;

				final_half_ls+=s_vector[i]*(f_1+f_2+f_3);
				
				

			}
		}
		__syncthreads();
		

     //3,nlevelsI,nlevelsF,sym%maxdegen,sym%maxdegen
		


	}
	if(valid)
	{
		half_ls[irootF] = final_half_ls;
	}	
	


}


///K-fast---method


__global__ void device_compute_1st_half_ls_flipped_dipole_shared_blocks_fast(const int dimenF,const int dimenI,const int j0dimen,const int jI,const int jF,
const int kF,const int* tauI_,const int* tauF_, const int* icorrI_,const int* icorrF_,const int startFblock,const int endFblock,const int ncontrF,
const double* __restrict__ dipole_me,
const double* vector,const double* __restrict__  threej,double*  half_ls)
{

	//const int irootF = blockIdx.x * blockDim.x + threadIdx.x;
	//double sq2 = 1.0/sqrt(2.0);
	volatile __shared__ double s_ls_factor[DIPOLE_BLOCK_SIZE];
	volatile __shared__ int s_icontrI[DIPOLE_BLOCK_SIZE];
	volatile __shared__ int s_tauI[DIPOLE_BLOCK_SIZE];
	volatile __shared__ double s_sigmaI[DIPOLE_BLOCK_SIZE];;
	int t_id = threadIdx.x;
	//int b_size = 32;//BLOCK_SIZE;
	int b_start = (threadIdx.x/WARP_SIZE)*WARP_SIZE;
	int w_id = threadIdx.x % WARP_SIZE;
	int icontrF, kI, tauI,tauF,sigmaF, sigmaI, dipole_idx;
	//These are o remove if statements
	int kI_kF_eq,tauF_tauI_neq,kI_kF_zero;
	double ls = 0.0,f3j=0.0,final_half_ls,sq_factor;
	int startF = (2*(kF-1)+1)*j0dimen;
	startF <0? startF=0: startF=startF;
	s_ls_factor[t_id]=0.0;
	int valid = 0;
	//__syncthreads();
	//for(int irootF=blockIdx.x * blockDim.x + threadIdx.x + startF; irootF < endF; irootF+=blockDim.x*gridDim.x){
	int irootF=blockIdx.x*blockDim.x + threadIdx.x + startF;
	//kF=-10000;
	icontrF=-10;
	tauF = 0;

	//for(irootF=blockIdx.x * blockDim.x + threadIdx.x + startF; irootF < endF; irootF+=blockDim.x*gridDim.x){	
		final_half_ls = 0.0;
	if(irootF<dimenF){
		//If we are out of range the we always acces the zeroth element
		icontrF =icorrF_[irootF];
		tauF  =  tauF_[irootF] & 1;
	}

	valid = ((icontrF >=startFblock) && (icontrF < endFblock));
	icontrF-=startFblock;

	sigmaF = (kF % 3)*tauF;
	kI = max(kF-1,0);
	
	for(kI=kI; kI<kF+2 && kI<(jI+1); kI++){
		int startI = min((2*(kI-1)+1)*j0dimen,0);
		int endI = (2*kI+1)*j0dimen;
		kI*kF!=0 ? sq_factor=SQRT2 : sq_factor=1.0;
		f3j  =  threej[jI + kI*(c_jmax+1) + (jF - jI + 1)*(c_jmax+1)*(c_jmax+1) + (kF - kI +1)*(c_jmax+1)*(c_jmax+1)*3];
		for(int b_irootI=startI; b_irootI < endI; b_irootI+=WARP_SIZE){

			int irootI = b_irootI+w_id;

			s_ls_factor[t_id]=0.0;

			if(irootI < b_irootI+endI){
				//Cache into memory
				s_tauI[t_id] = tauI_[irootI] & 1;
				s_icontrI[t_id] = icorrI_[irootI];
				
				sigmaI = (kI % 3)*tauI;
				sigmaI = 2*(~(sigmaI+kI) & 1)-1;


				s_sigmaI[t_id] = double(sigmaI);
				//Perform calculations
				kI_kF_eq = (kF==kI); 	  
				kI_kF_zero = ((kI*kF) != 0); 

				ls =f3j*vector[irootI]*sq_factor;
				s_ls_factor[t_id] = ls;//1.0;//ls;
			}
					//__threadfence();
				
					//__syncthreads();
			for(int i = 0; i < WARP_SIZE; i++){

				ls= s_ls_factor[i + b_start];

				if(ls==0.0) continue;
				tauI=s_tauI[i + b_start];
				tauF_tauI_neq = (tauF!=tauI);

				dipole_idx=2*(kI_kF_eq)+ (!kI_kF_eq)*(!tauF_tauI_neq);

				ls *= (tauF-tauI)*(kI_kF_eq) + (tauF-tauI)*(kF-kI)*(!kI_kF_eq)*( tauF_tauI_neq) - (!kI_kF_eq)*(!tauF_tauI_neq);

				if(valid) final_half_ls+= s_sigmaI[i + b_start]*ls*dipole_me[icontrF + s_icontrI[i + b_start]*ncontrF + dipole_idx*dip_stride_1*ncontrF];
					
			} 
	
	
			//s_ls_factor[t_id]=0.0;
							
		}
	}
	

	if(irootF<dimenF && valid){
		final_half_ls *= double(2*(~(sigmaF) & 1)-1);
		half_ls[irootF] = final_half_ls;
	}
}



__global__ void device_sort_vector(const int dimenI,const int* N_, const double* c_vec_,double* vec_){

	int irootI = blockIdx.x * blockDim.x + threadIdx.x; 

	if(irootI < dimenI){
		vec_[irootI] = c_vec_[N_[irootI]];
		
	}	
}


__global__ void device_unsort_vector(const int dimenI,const int* N_, const double* c_vec_,double* vec_){

	int irootI = blockIdx.x * blockDim.x + threadIdx.x; 

	if(irootI < dimenI){
		vec_[N_[irootI]] = c_vec_[irootI];
		
	}	
}


//This relates to TROVEs method of expanding vectors to the primitive basis set
__global__ void device_expand_vectors(const int dimenI,const int igammaI,const int maxcoeff,const int idegI,const int sdeg,const int Ntot,const int* ijterms_, const int* icontr_,const int* N_,const double* repres_, const double* vecI_,double* vec_)
{
	int irootI = blockIdx.x * blockDim.x + threadIdx.x; 
	if(irootI < dimenI)
	{
	
		int irow;
		int iterm;
		int nelem;
		int isrootI;
		int ib; 
		double dtemp0 = 0.0;
		irow = icontr_[irootI];

		ib = icontr_[irootI + dimenI];
	
		iterm = ijterms_[irow + igammaI*maxcoeff];
	
		nelem = N_[irow];
		
		for(int i = 0; i < nelem; i++)
		{
			isrootI = iterm+i;
			dtemp0 +=  vecI_[isrootI]*repres_[isrootI + idegI*Ntot + ib*sdeg*Ntot];
		}
	
		vec_[irootI] = dtemp0;
	}


}



void cleanup_gpu(){
	

}

