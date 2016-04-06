#include "cublas_v2.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cmath>
#include <cstdio>
#include <vector>
//Kepler definitions
#ifdef KEPLER
	#define CORRELATE_BLOCK_SIZE 512
	#define DIPOLE_BLOCK_SIZE 128
	#define MAX_STREAMS 32
	#define MAX_VECTORS 64
#else
	#define CORRELATE_BLOCK_SIZE 512
	#define DIPOLE_BLOCK_SIZE 128
	#define MAX_STREAMS 16
	#define MAX_VECTORS 64
#endif

//Math constants

//#define DEBUG
#define WARP_SIZE 32
#define PI 3.14159265359
#define AVOGNO 6.0221415e+23
#define PLANCK 6.62606896e-27
#define VELLGT 2.99792458e+10
#define BOLTZ 1.380658e-16
#define ACOEF 3.13618923e-7
#define SQRT2 0.70710678118

/*---------Basis set information---------------*/
std::vector<int> J;
std::vector<int> dimen;
std::vector< std::vector<int> > h_k_blocks;
std::vector< std::vector<int> > k_start;
std::vector< std::vector<int > > nlevelI;
std::vector< std::vector<int > > nlevelF;
std::vector<int> max_sym_coeff;
class GainGpuData{
public:
	std::vector<int*> K;
	std::vector<int*> icorr;
	std::vector<int*> tau;



	std::vector<int*> k_block_size;

	std::vector<int*> icontr;
	std::vector< std::vector<int*> > irr_N;
	std::vector< std::vector<double*> > irr_repres;
	std::vector<int*> ijterms;


	std::vector<cudaStream_t> ddot_stream;
	std::vector<cudaStream_t> half_ls_stream;

	double* corr_vector;
	std::vector< double* > vector_f;
	std::vector< double* > corr_vector_f;
	std::vector< double* > ls;
	std::vector< std::vector< double*> >  g_half_ls;
	std::vector< std::vector< double* > > wigner;

	double* g_threeJ;
	double* cuda_dipole_me;
	size_t gpu_free_space;
	std::vector< cublasHandle_t > handle;
	cublasStatus_t stat;
};

std::vector<GainGpuData> workers;


int cuda_dimenmax;
int nsizemax;






int nJ;


int degen_max;
int sym_nrepres;

bool spread_dipoles = false;
bool big_dipole = false;


//Dot product streams

int num_procs;
int num_gpus_per_proc=-1;

int verbose;

//Classes
struct DipolePiece{
	double* dipole_me;
	int startF;
	int endF;
	int ncontrF;
	size_t size;
};


//Reads the dipole_ flipped
class GainDipole{
private:
	std::vector<DipolePiece> dipole_me;
	size_t dipole_size;
	int MaxContracts;
	int num_blocks;
	size_t max_size;
	size_t max_memory;
	int m_parts;
public:		
	GainDipole(size_t pmax_memory, int parts) : max_memory(pmax_memory){};

	void SplitDipole(double * dipole_me,int ncontr_t);
	//const double* GetDipole(){ return dipole_me;}
	size_t GetDipoleSize(){ return dipole_size;}
	size_t GetDipoleSizeBytes(){return dipole_size*sizeof(double);}
	int GetNumBlocks(){return num_blocks;};

	const double* GetDipolePiece(int block){ return dipole_me.at(block).dipole_me;};
	int GetDipoleStart(int block){ return dipole_me.at(block).startF;};
	int GetDipoleEnd(int block){ return dipole_me.at(block).endF;};
	int GetDipoleNcontr(int block){ return dipole_me.at(block).ncontrF;};
	size_t GetDipoleSize(int block){ return dipole_me.at(block).size;};
	size_t GetMaxDipoleSize(){ size_t dip_size=0; for(int i = 0; i < dipole_me.size(); i++){ max_size = max((int)dipole_me.at(i).size,(int)max_size);} return max_size;}
	size_t GetMaxContracts(){ return MaxContracts;};
	void Cleanup();
	~GainDipole();

};


void GainDipole::SplitDipole(double * pdipole_me,int ncontr_t){
	char buff[20];
	int imu,imu_t;
	size_t matsize,rootsize,rootsize_t,rootsize2;

	rootsize2 = ncontr_t*ncontr_t;
	matsize = rootsize2*3;
	dipole_size = matsize;
	size_t total_mat_bytes = matsize*sizeof(double);
	if( m_parts == -1)
		num_blocks = (total_mat_bytes/max_size) + 1;
	else
		num_blocks = m_parts;

	rootsize  = ncontr_t*(ncontr_t+1)/num_blocks;
	int n_contr_block = ceil(float(ncontr_t)/float(num_blocks));
	int cur_block_size = 0;
	int startF=0,endF=0,ncontrF=0;
	//Allocate X parts
	//dipole_block.dip_block=new FDipole_block[parts];
	//dipole_block.parts=parts;	
	printf("Flipping dipole...and blocking into %i pieces\n",num_blocks);
	for(int blocks = 0; blocks < num_blocks; blocks++){
		startF=blocks*n_contr_block;
		
		endF = (blocks+1)*n_contr_block;
		endF = std::min(endF,ncontr_t);
		ncontrF = endF-startF;
		dipole_me.push_back(DipolePiece());
		//Allocate the dipole
		dipole_me.back().startF = startF;
		dipole_me.back().endF = endF;
		dipole_me.back().dipole_me = new double[ncontrF*ncontr_t*3];
		dipole_me.back().size = sizeof(double)*size_t(ncontrF)*size_t(ncontr_t)*3l;
		printf("Size of block %i is %zu\n",blocks,dipole_me.back().size);
		dipole_me.back().ncontrF = ncontrF;
		printf("startF =%i, block number = %i, n_contr_block = %i\n",startF,blocks,ncontrF);
		for(int i = 0; i < ncontr_t; i++)
			for(int f = startF; f < endF; f++)
				for(int k = 0; k < 3; k++)
					dipole_me.back().dipole_me[f-startF + i*ncontrF + k*ncontr_t*ncontrF] = pdipole_me[f + i*ncontr_t + k*ncontr_t*ncontr_t]; 

	}


	//delete[] temp_dipole;
	printf("done!\n");
}

void GainDipole::Cleanup(){
	for(int i = 0; i < dipole_me.size(); i++){
		if(dipole_me.at(i).dipole_me != NULL)
			delete [] dipole_me.at(i).dipole_me;
		dipole_me.at(i).dipole_me = NULL;
	}
}


GainDipole::~GainDipole(){
//	if(dipole!=NULL)
//		delete[] dipole;

}




GainDipole* gpu_dipole;





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


__global__ void device_expand_vectors(const int dimenI,const int igammaI,const int maxcoeff,const int idegI,const int sdeg,const int Ntot,const int* ijterms_, const int* icontr_,const int* N_,const double* repres_, const double* vecI_,double* vec_);


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
const double* vector,const double*  threej,double*  half_ls);


void CheckCudaError(const char* tag){
  // check for error
  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
  {
    // print the CUDA error message and exit
    printf("[%s] CUDA error: %s\n", tag,cudaGetErrorString(error));
    cudaDeviceReset();

    exit(-1);
  }
};


extern "C" void gpu_pin_vector_memory_(double* vector,int* n){
	if(cudaSuccess!= cudaHostRegister(vector,sizeof(double)*size_t(*n),cudaHostAllocPortable)){
		printf("Problem pinning vector memory!!\n");
		exit(0);
	}
}

extern "C" void gpu_unpin_vector_memory_(double* vector){
	if(cudaSuccess!= cudaHostUnregister(vector)){
		printf("Problem unpinning vector memory!!\n");
		exit(0);
	}
}





void initialise_gpu_(int gpu_id,int* dimenmax,int * verbose_,int* jmax_,int* nJ_,int*degen_max_,int* maxsym_,int* num_procs_){
		cudaSetDevice(gpu_id);

		workers.push_back(GainGpuData());		
		
		printf("%i degen=%i nJ=%i jmax=%i,sym_nrepres=%i nsizemax=%i\n",gpu_id,degen_max,*nJ_,*jmax_,*maxsym_,nsizemax);

		cudaFree(0);
			//Get device properties
		cudaDeviceProp devProp;
		cudaGetDeviceProperties(&devProp, 0);
			//Now we compute total bytes
		workers.back().gpu_free_space = size_t(double(devProp.totalGlobalMem)*0.95);

		//Degen max
		


		if(verbose > 4 ) printf("We have %zu bytes to play with\n",workers.back().gpu_free_space);	
	
		if(verbose > 3) printf("\nInitializing vectors in GPU memory (dimen max = %i).....",cuda_dimenmax);		

		if(cudaSuccess != cudaMalloc((void**)&workers.back().corr_vector,sizeof(double)*cuda_dimenmax)){
			printf("Error allocating corralated vectors!\n");
			exit(0);	
		}
		printf(".We have %i streams to play with...",num_procs);
		for(int i = 0; i < num_procs; i++){


			workers.back().ddot_stream.push_back(NULL);
			cudaStreamCreate(&workers.back().ddot_stream.back());
			workers.back().corr_vector_f.push_back(NULL);
			workers.back().vector_f.push_back(NULL);
			workers.back().ls.push_back(NULL);
			if(cudaSuccess != cudaMalloc((void**)&workers.back().corr_vector_f.back(),sizeof(double)*cuda_dimenmax)){
				printf("Error allocating corralated vectors!\n");
				exit(0);	
			}
			if(cudaSuccess != cudaMalloc((void**)&workers.back().vector_f.back(),sizeof(double)*nsizemax)){
				printf("Error allocating nocorrelates vectors!\n");
				exit(0);	
			}
			if(cudaSuccess != cudaMalloc((void**)&workers.back().ls.back(),sizeof(double)*degen_max*degen_max)){
				printf("Error allocating corralated vectors!\n");
				exit(0);	
			}
			workers.back().gpu_free_space -= sizeof(double)*nsizemax;
			workers.back().gpu_free_space -= sizeof(double)*cuda_dimenmax;
			workers.back().gpu_free_space -= sizeof(double)*degen_max;
		}

		for(int i = 0; i < nJ; i++){
			workers.back().g_half_ls.push_back(std::vector<double*>());
			for(int j = 0; j < degen_max; j++){
				workers.back().g_half_ls.back().push_back(NULL);
				if(cudaSuccess != cudaMalloc((void**)&workers.back().g_half_ls.back().back(),sizeof(double)*cuda_dimenmax)){
					printf("Error allocating half_ls vector!\n");
					exit(0);	
				}				
				workers.back().gpu_free_space -= sizeof(double)*cuda_dimenmax;
			}
		}


		printf("Initialized vectors\n");

		printf("Initializing cuBlas\n");
		for(int i = 0; i < num_procs; i++){
			workers.back().handle.push_back(0);
			workers.back().stat = cublasCreate(&workers.back().handle.back());
			if (workers.back().stat != CUBLAS_STATUS_SUCCESS) {
				printf ("CUBLAS initialization failed\n");
				return;
			}
			cublasSetPointerMode(workers.back().handle[i], CUBLAS_POINTER_MODE_DEVICE);
			cublasSetStream(workers.back().handle[i],workers.back().ddot_stream[i]);
		}
		
		for(int i = 0; i < MAX_STREAMS; i++){

			workers.back().half_ls_stream.push_back(NULL);
			cudaStreamCreate(&workers.back().half_ls_stream.back());
		}


		workers.back().gpu_free_space -= sizeof(double)*cuda_dimenmax;

		host_jmax=*jmax_;
		CheckCudaError("Alloc vectors");		
		cudaDeviceSynchronize();







		cudaMemcpyToSymbol (g_maxdegen,&degen_max , sizeof(int) );		
		CheckCudaError("maxdegen");

		cudaMemcpyToSymbol (c_jmax,&host_jmax , sizeof(int) );		
		CheckCudaError("Jmax symbol");
		cudaMemcpyToSymbol ( g_symnrepres,&nsizemax, sizeof(int) );	
		CheckCudaError("Nrepres");
		if(verbose>3) printf("done!\n");
}


extern "C" void initialize_gain_(int* dimenmax,int * verbose_,int* jmax_,int* nJ_,int*degen_max_,int* maxsym_,int* num_procs_){
	char* pNumGPU;
	verbose = *verbose_;
	cuda_dimenmax = *dimenmax;
	nsizemax = *dimenmax;
	degen_max = *degen_max_;

	nJ=*nJ_;
	num_procs=*num_procs_;
	sym_nrepres = *maxsym_;


	if(pNumGPU != NULL){
		num_gpus_per_proc = atoi(pNumGPU);
	}else{
		cudaGetDeviceCount(&num_gpus_per_proc);
		
	}
	printf("We are working with %d GPUs\n",num_gpus_per_proc);
	for(int i = 0; i < num_gpus_per_proc; i++){
		initialise_gpu_(i,dimenmax,verbose_,jmax_,nJ_,degen_max_,maxsym_,num_procs_);
	}
}

extern "C" void transfer_dipole_(double* dipole_me,int* ncontr_t_){
		cudaSetDevice(0);
		int ncontr_t = *ncontr_t_;
		int max_elem = ncontr_t*ncontr_t*3;
		int ncontr_t_2 = ncontr_t*ncontr_t;
		size_t max_mgpu = 0;

		for(int i = 0; i < workers.size(); i++){
			max_mgpu += workers[i].gpu_free_space;
		}

		if (sizeof(double)*size_t(max_elem) > workers[0].gpu_free_space){
			printf(" Dipole too big, splitting!\n");
			gpu_dipole = new GainDipole(workers[0].gpu_free_space,-1);
			//exit(0);
		}
		
		if(sizeof(double)*size_t(max_elem) > max_mgpu){
			printf("Dipole too big for all GPUs\n");
			gpu_dipole = new GainDipole(workers[0].gpu_free_space,-1);
			big_dipole = true;
			//exit(0);
		}else{
			gpu_dipole = new GainDipole(workers[0].gpu_free_space,workers.size());
			spread_dipoles = true;
		} 

		printf("Transfering dipoles");
		gpu_dipole->SplitDipole(dipole_me,*ncontr_t_);

		for(int i = 0; i < workers.size(); i++){
			cudaSetDevice(i);
			if(cudaSuccess != cudaMalloc((void**)&workers.at(i).cuda_dipole_me,gpu_dipole->GetMaxDipoleSize())){
				printf("Error dipole!\n");
				exit(0);			
			}
			//if(spread_dipole){	
				cudaMemcpy(workers.at(i).cuda_dipole_me,gpu_dipole->GetDipolePiece(i),gpu_dipole->GetDipoleSize(i),cudaMemcpyHostToDevice);			
			//}
			if(!spread_dipoles)
				break;
			
		}
		
		
		if(!big_dipole)
			gpu_dipole->Cleanup();

		cudaDeviceSynchronize();
		CheckCudaError("Dipole");

		ncontr_t_2 = ncontr_t*ncontr_t;

		for(int i = 0; i < workers.size(); i++){
			cudaSetDevice(i);
			cudaMemcpyToSymbol (dip_stride_1,&ncontr_t , sizeof(int) );
			cudaDeviceSynchronize();
			CheckCudaError("Symbol");
			cudaMemcpyToSymbol (dip_stride_2,&ncontr_t_2 , sizeof(int) );
		}



/*
		if(cudaSuccess != cudaMalloc((void**)&cuda_dipole_me,sizeof(double)*size_t(max_elem))){
			printf("Error dipole!\n");
			exit(0);			
		}
		if(verbose>3) printf("\nTransfering Dipole to GPU.....");
		double* temp_dipole= new double[ncontr_t*ncontr_t*3];
		for(int i = 0; i < ncontr_t; i++)
			for(int f = 0; f < ncontr_t; f++)
				for(int k = 0; k < 3; k++){
					temp_dipole[f + i*ncontr_t + k*ncontr_t*ncontr_t] = dipole_me[i + f*ncontr_t + k*ncontr_t*ncontr_t]; 
				}
		
		cudaMemcpy(cuda_dipole_me,temp_dipole,sizeof(double)*size_t(max_elem),cudaMemcpyHostToDevice);

		cudaDeviceSynchronize();
		CheckCudaError("Dipole");
		for(int i = 0; i < workers.size(); i++){
			cudaMemcpyToSymbol (dip_stride_1,&ncontr_t , sizeof(int) );
			cudaDeviceSynchronize();
			CheckCudaError("Symbol");
			cudaMemcpyToSymbol (dip_stride_2,&ncontr_t_2 , sizeof(int) );
		}
		if(verbose>3) printf("done!\n");
		delete[] temp_dipole;
		gpu_free_space -= sizeof(double)*size_t(max_elem);
		
*/



}

extern "C" void transfer_threej_(double* three_J){
	if(verbose>3) printf("\nTransfering threej to GPUs.....");
	for(int i = 0; i < workers.size(); i++){
		cudaSetDevice(i);	
		cudaMalloc((void**)&workers.at(i).g_threeJ,size_t((host_jmax+1)*(host_jmax+1)*3*3)*sizeof(double));
		workers.at(i).gpu_free_space -= (host_jmax+1)*(host_jmax+1)*3*3*sizeof(double);
		cudaMemcpy(workers.at(i).g_threeJ,three_J,size_t((host_jmax+1)*(host_jmax+1)*3*3)*sizeof(double),cudaMemcpyHostToDevice);
		if(verbose>3) printf("done!\n");
		CheckCudaError("ThreeJ");
	}	
}

__host__ int* transfer_basis_parts(const char* name,int* item,size_t size){
	int* temp_ptr;
	if(cudaSuccess != cudaMalloc((void**)&temp_ptr,size)){
			printf("Error allocating for %s!\n",name);
			exit(0);	
	}
	cudaMemcpy(temp_ptr,item,size,cudaMemcpyHostToDevice);	
	cudaDeviceSynchronize();
	//gpu_free_space -= size;
	return temp_ptr;	
	
}

__host__ double* transfer_basis_parts_d(const char* name,double* item,size_t size){
	double* temp_ptr;
	if(cudaSuccess != cudaMalloc((void**)&temp_ptr,size)){
			printf("Error allocating for %s!\n",name);
			exit(0);	
	}
	cudaMemcpy(temp_ptr,item,size,cudaMemcpyHostToDevice);	
	cudaDeviceSynchronize();
	CheckCudaError(name);
	//gpu_free_space -= size;
	return temp_ptr;	
	
}

extern "C" double gpu_dot_(int*proc_,int* n_,int* indF,int* idegI,double* vector){
	/*
	int proc=*proc_;
	double result;
	int n = *n_;
	//for(int i = 0; i < n; i++)
	//	printf("%16.4E  %16.4E\n",x[i],y[i]);
	//cublasSetVector(n,sizeof(double),vector,1,corr_vector_f,1);	
	cudaMemcpyAsync(corr_vector_f[proc],vector,sizeof(double)*size_t(n), cudaMemcpyHostToDevice,ddot_stream[proc]);
	
	cublasDdot (handle[proc],n,g_half_ls[*indF -1][*idegI -1],1,corr_vector_f[proc],1,ls[proc]);
	cudaMemcpyAsync(&result,ls[proc],sizeof(double), cudaMemcpyDeviceToHost,ddot_stream[proc]);
	cudaStreamSynchronize(ddot_stream[proc]);
	*/

	return -1.0;//result;

}

extern "C" void correlate_vectors(int*proc_,int* n_,int* indF_,int* igammaF_,int* idegI_,int* idegF_,int* sdeg_,int* Ntot_,double* vector){
	
	int proc=*proc_;
	int n = *n_;
	int indF=*indF_-1;
	int idegI=*idegI_-1;
	int idegF=*idegF_-1;
	int dimenF = dimen[indF];
	int sym_coeff = max_sym_coeff[indF];	
	int igammaF=*igammaF_-1;
	int sdeg=*sdeg_;
	int Ntot=*Ntot_;
	int choose_gpu = proc % num_gpus_per_proc;
	cudaMemcpyAsync(workers.at(choose_gpu).vector_f[proc],vector,sizeof(double)*size_t(n), cudaMemcpyHostToDevice,workers.at(choose_gpu).ddot_stream[proc]);
	//cudaDeviceSynchronize();
	//CheckCudaError("MemcpyAsync");
//__global__ void device_expand_vectors(
		//const int dimenI,
		//const int igammaI,
		//const int maxcoeff,
		//const int idegI,
		//const int sdeg,
		//const int Ntot,
		//const int* ijterms_,
		// const int* icontr_,
		//const int* N_,
		//const double* repres_, 
		//const double* vecI_,
		//double* vec_);
	int gridSize = (int)ceil((float)dimenF/CORRELATE_BLOCK_SIZE);
	device_expand_vectors<<<gridSize,CORRELATE_BLOCK_SIZE,0,workers.at(choose_gpu).ddot_stream[proc]>>>(
		dimenF,
		igammaF,
		sym_coeff,
		idegI,
		sdeg,
		Ntot,
		workers.at(choose_gpu).ijterms[indF],
		workers.at(choose_gpu).icontr[indF],
		workers.at(choose_gpu).irr_N[indF][igammaF],
		workers.at(choose_gpu).irr_repres[indF][igammaF],
		workers.at(choose_gpu).vector_f[proc],
		workers.at(choose_gpu).corr_vector_f[proc]);	
}


extern "C" void corr_gpu_dot_(int*proc_,int* n_,int* indF_,int* igammaF_,int* idegI_,int* idegF_,int* sdeg_,int* Ntot_,double* vector){
	
	int proc=*proc_;
	int n = *n_;
	int indF=*indF_-1;
	int idegI=*idegI_-1;
	int idegF=*idegF_-1;
	int dimenF = dimen[indF];
	int sym_coeff = max_sym_coeff[indF];	
	int igammaF=*igammaF_-1;
	int sdeg=*sdeg_;
	int Ntot=*Ntot_;
	
	int choose_gpu = proc % num_gpus_per_proc;
	cudaSetDevice(choose_gpu);
	//printf("proc=%i idegI=%i idegF = %i sdeg=%i Ntot=%i\n",proc,idegI,idegF,sdeg,Ntot);
		
	//for(int i = 0; i < n; i++)
	//	printf("%16.4E  %16.4E\n",x[i],y[i]);
	//cublasSetVector(n,sizeof(double),vector,1,corr_vector_f,1);	
	cudaMemcpyAsync(workers.at(choose_gpu).vector_f[proc],vector,sizeof(double)*size_t(n), cudaMemcpyHostToDevice,workers.at(choose_gpu).ddot_stream[proc]);
	//cudaDeviceSynchronize();
	//CheckCudaError("MemcpyAsync");
//__global__ void device_expand_vectors(
		//const int dimenI,
		//const int igammaI,
		//const int maxcoeff,
		//const int idegI,
		//const int sdeg,
		//const int Ntot,
		//const int* ijterms_,
		// const int* icontr_,
		//const int* N_,
		//const double* repres_, 
		//const double* vecI_,
		//double* vec_);
	int gridSize = (int)ceil((float)dimenF/CORRELATE_BLOCK_SIZE);
	device_expand_vectors<<<gridSize,CORRELATE_BLOCK_SIZE,0,workers.at(choose_gpu).ddot_stream[proc]>>>(
		dimenF,
		igammaF,
		sym_coeff,
		idegI,
		sdeg,
		Ntot,
		workers.at(choose_gpu).ijterms[indF],
		workers.at(choose_gpu).icontr[indF],
		workers.at(choose_gpu).irr_N[indF][igammaF],
		workers.at(choose_gpu).irr_repres[indF][igammaF],
		workers.at(choose_gpu).vector_f[proc],
		workers.at(choose_gpu).corr_vector_f[proc]);	

	//cudaDeviceSynchronize();
	//CheckCudaError("expand");
	cublasDdot (workers.at(choose_gpu).handle[proc],dimenF,workers.at(choose_gpu).g_half_ls[indF ][idegI ],1,workers.at(choose_gpu).corr_vector_f[proc],1,workers.at(choose_gpu).ls[proc]+ (idegI + idegF*degen_max));
	//cudaDeviceSynchronize();
	//CheckCudaError("Ddot");
	//cudaMemcpyAsync(&result,ls[proc],sizeof(double), cudaMemcpyDeviceToHost,ddot_stream[proc]);
	//cudaStreamSynchronize(ddot_stream[proc]);
	
	//return result;

}

extern "C" void retrieve_linestrength_(int* proc,double* ls_){
	int choose_gpu = *proc % num_gpus_per_proc;
	cudaSetDevice(choose_gpu);
	cudaMemcpyAsync(ls_,workers.at(choose_gpu).ls[*proc],sizeof(double)*degen_max*degen_max, cudaMemcpyDeviceToHost,workers.at(choose_gpu).ddot_stream[*proc]);
	//cudaDeviceSynchronize();
	//CheckCudaError("Retrieve");
	cudaStreamSynchronize(workers.at(choose_gpu).ddot_stream[*proc]);
}


extern "C" size_t get_memory_allocation(){}


extern "C" void transfer_basis_set_(int* j_, int* k_, int * ktau_,int* icorr_,int* icontr_, int* dimen_,int* ijterm_,int* maxsym_,int* do_rotsym){

	//push_back_J	
	int jval=*j_;
	J.push_back(jval);
	int maxcontracts = *dimen_;

	dimen.push_back(maxcontracts);
	
	for(int i = 0; i < workers.size(); i++){
		cudaSetDevice(i);
		max_sym_coeff.clear();
		h_k_blocks.clear();
		k_start.clear();
		if(verbose>3) printf("\nTransfering basis-set for J=%d to GPU %i....\n.",*j_,i);
		//Transfer K
		workers.at(i).K.push_back(transfer_basis_parts("K",k_,sizeof(int)*maxcontracts));
		//if(verbose>4) printf("[%d] transfering K",jval);
		//Transfer correlation
		workers.at(i).icorr.push_back(transfer_basis_parts("icorr",icorr_,sizeof(int)*maxcontracts));
		//if(verbose>4) printf("[%d] transfering icorr",jval);
		//Transfer Tau	
		if(*do_rotsym == 0){
			int* t_tau = new int[maxcontracts];
			for(int k = 0; k < maxcontracts; k++){
				t_tau[k] = ktau_[k] % 2;
				//if(verbose>4) printf("tau[%i]=%i\n",i,t_tau[i]);
			}
			workers.at(i).tau.push_back(transfer_basis_parts("tau",t_tau,sizeof(int)*maxcontracts));
			delete[] t_tau;
		}else{
			printf("Normal ktau\n");
			workers.at(i).tau.push_back(transfer_basis_parts("tau",ktau_,sizeof(int)*maxcontracts));
		}

		//Perform k block analysis
		if(*do_rotsym == 0){
			//Now we need to perform some K anaqlysis
			h_k_blocks.push_back(std::vector<int>());	
			k_start.push_back(std::vector<int>());
	
			//printf("doing K block stuff\n");
			int* temp_k_block_size = new int[jval+1];
			//// K bloack analysis

			int last_k_b = 0;
			int ksize=0;

			k_start.back().push_back(0);
			for(int k = 0; k < maxcontracts; k++){
				int t_k = k_[k];
				//printf("K = %i last_k_b=%i ksize=%i \n",t_k,last_k_b,ksize);
				if(last_k_b != t_k){
					temp_k_block_size[last_k_b]=ksize;
					h_k_blocks.back().push_back(ksize);
					k_start.back().push_back(k_start.back()[last_k_b]+ksize);
					ksize=0;
					last_k_b=t_k;
				}

				ksize++;

			}
			h_k_blocks.back().push_back(ksize);
			temp_k_block_size[jval]=ksize;
			if(verbose>3) printf("\n\n---------K block analysis for J=%d------------\n",jval);
			for(int k =0; k < jval+1; k++){
				if(verbose>3) printf("K block %i start: %i end %d size: %i \n",k,k_start.back()[k],k_start.back()[k]+h_k_blocks.back()[k]-1,h_k_blocks.back()[k]);
				fflush(0);
			}
			printf("\n");
			//exit(0);
			cudaDeviceSynchronize();
			CheckCudaError("Basis-set");

			workers.at(i).k_block_size.push_back(transfer_basis_parts("kblcok",temp_k_block_size,sizeof(int)*(jval+1)));
			delete[] temp_k_block_size;
		}
		//icontr2icase
		workers.at(i).icontr.push_back(transfer_basis_parts("icontr",icontr_,sizeof(int)*(maxcontracts)*(2)));

		//ijterm
		workers.at(i).ijterms.push_back(transfer_basis_parts("ijterm",ijterm_,sizeof(int)*(sym_nrepres)*(*maxsym_)));
	
		workers.at(i).irr_N.push_back(std::vector<int*>());
		workers.at(i).irr_repres.push_back(std::vector<double*>());
	
		max_sym_coeff.push_back(*maxsym_);
	
		workers.at(i).wigner.push_back(std::vector<double*>());


		if(!spread_dipoles)
			break;
	}
	nlevelI.push_back(std::vector<int>());
	nlevelF.push_back(std::vector<int>());
	if(verbose>3) printf("done!\n");
	
}

extern "C" void transfer_wigner_(double* wigner_,int* nlevelI_,int* nlevelF_){
	printf("nlevelI = %i nlevelF = %i degen_max = %i\n",*nlevelI_,*nlevelF_,degen_max);
	int nlevelI__ = *nlevelI_;
	int nlevelF__ = *nlevelF_;
	nlevelI.back().push_back(nlevelI__);
	nlevelF.back().push_back(nlevelF__);
	for(int i = 0; i < workers.size(); i++){
		cudaSetDevice(i);
		workers.at(i).wigner.back().push_back(transfer_basis_parts_d("wigner",wigner_,sizeof(double)*3*nlevelI__*nlevelF__*degen_max*degen_max));
	}	
	printf("%p\n",workers[0].wigner.back().back());
}


extern "C" void transfer_inflation_(int* N, double* repres,int* maxsym_,int* Ntot_,int* degen_,int* mat_size_){
	printf("Ntot = %i, degen=%i mat_size=%i max_sym_=%i\n",*Ntot_,*degen_,*mat_size_,*maxsym_);
	int Ntot=*Ntot_;
	int degen = *degen_;
	int mat_size=*mat_size_;
	int max_sym=*maxsym_;
	for(int i = 0; i < workers.size(); i++){
		cudaSetDevice(i);
		workers.at(i).irr_N.back().push_back(transfer_basis_parts("N",N,sizeof(int)*max_sym));
	//for(int i = 0; i < max_sym; i++)
	//	printf("N[%i]=%i\n",i,N[i]);
	//exit(0);
		workers.at(i).irr_repres.back().push_back(transfer_basis_parts_d("repres",repres,sizeof(double)*Ntot*degen*mat_size));
	}
	
}


void compute_gpu_half_linestrength_small(int* indI_,int* indF_,int* ideg_,int* n_,double* vector,double* half_ls){

			int indI = *indI_ -1;
			int indF = *indF_ -1;
			int ideg = *ideg_ -1;
			//memcpy(
			//	cudaMemcpyAsync(workers.at(0).vector_f[0],vector,sizeof(double)*size_t(n), cudaMemcpyHostToDevice,workers.at(0).ddot_stream[proc]);
				//cudaMemcpyAsync(workers.at(choose_gpu).vector_f[proc],vector,sizeof(double)*size_t(n), cudaMemcpyHostToDevice,workers.at(choose_gpu).ddot_stream[proc]);
			//cudaMemcpy(workers[0].vector_f,vector,sizeof(double)*cuda_dimenmax,cudaMemcpyHostToDevice);
			

			

			cudaDeviceSynchronize();
			CheckCudaError("Memcpy");
			int jF = J[indF];
			int jI = J[indI];
			int stream_id;
			//printf("jI %i jF %i\n",jI,jF);
			int total_k_blocks = jF + 1;
			for(int k = 0; k < total_k_blocks; k++){
				stream_id = k % MAX_STREAMS;
				//printf("K=%i startK = %i\n",k,k_start[indI][max(k-1,0)]);
                                fflush(0);

				int half_grid_size = (int)ceil((float)(h_k_blocks[indF][k])/(float)DIPOLE_BLOCK_SIZE);
				//printf("half_grid params: b:%i t:%i N:%i\n",half_grid_size,DIPOLE_BLOCK_SIZE,DIPOLE_BLOCK_SIZE*half_grid_size);
				device_compute_1st_half_ls_flipped_dipole_shared_nontrove<<<half_grid_size,DIPOLE_BLOCK_SIZE,0,workers[0].half_ls_stream[stream_id]>>>(h_k_blocks[indF][k],dimen[indI],jI,
											jF,k,workers[0].tau[indI],workers[0].tau[indF],workers[0].icorr[indI],workers[0].icorr[indF],
											k_start[indF][k],k_start[indI][max(k-1,0)],workers[0].k_block_size[indI],
											workers[0].cuda_dipole_me,
											workers[0].corr_vector,workers[0].g_threeJ,
											workers[0].g_half_ls[indF][ideg]);
				//cudaDeviceSynchronize();
				//CheckCudaError("K-block");
			}
			cudaDeviceSynchronize();
			CheckCudaError("Half_ls");
			//Copy to devices
			for(int i = 1; i < workers.size(); i++){
				cudaMemcpyPeer(workers[i].g_half_ls[indF][ideg],i,workers[0].g_half_ls[indF][ideg],0,sizeof(double)*cuda_dimenmax);
			}


			cudaMemcpy(half_ls,workers[0].g_half_ls[indF][ideg],sizeof(double)*cuda_dimenmax,cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();

}

void compute_gpu_half_linestrength_spread(int* indI_,int* indF_,int* ideg_,int* n_,double* vector,double* half_ls){

			int indI = *indI_ -1;
			int indF = *indF_ -1;
			int ideg = *ideg_ -1;
			//memcpy(


			for(int i = 0; i < workers.size(); i++){
				cudaSetDevice(i);
				cudaMemset(workers[i].g_half_ls[indF][ideg], 0,sizeof(double)*cuda_dimenmax);
					//cudaMemcpyAsync(workers.at(choose_gpu).vector_f[proc],vector,sizeof(double)*size_t(n), cudaMemcpyHostToDevice,workers.at(choose_gpu).ddot_stream[proc]);
				cudaMemcpy(workers[i].corr_vector,vector,sizeof(double)*cuda_dimenmax,cudaMemcpyHostToDevice);
				cudaDeviceSynchronize();
				CheckCudaError("Memcpy");
				int jF = J[indF];
				int jI = J[indI];
				int stream_id;
				//printf("jI %i jF %i\n",jI,jF);
				int total_k_blocks = jF + 1;
				for(int k = 0; k < total_k_blocks; k++){
					stream_id = k % MAX_STREAMS;
					//printf("K=%i startK = %i\n",k,k_start[indI][max(k-1,0)]);
		                        fflush(0);

					int half_grid_size = (int)ceil((float)(h_k_blocks[indF][k])/(float)DIPOLE_BLOCK_SIZE);
					//printf("half_grid params: b:%i t:%i N:%i\n",half_grid_size,DIPOLE_BLOCK_SIZE,DIPOLE_BLOCK_SIZE*half_grid_size);
					device_compute_1st_half_ls_flipped_dipole_shared_block_nontrove<<<half_grid_size,DIPOLE_BLOCK_SIZE,0,workers[i].half_ls_stream[stream_id]>>>(h_k_blocks[indF][k],dimen[indI],jI,
												jF,k,workers[i].tau[indI],workers[i].tau[indF],workers[i].icorr[indI],workers[i].icorr[indF],
												k_start[indF][k],k_start[indI][max(k-1,0)],
												gpu_dipole->GetDipoleStart(i),
												gpu_dipole->GetDipoleEnd(i),
												gpu_dipole->GetDipoleNcontr(i),
												workers[i].k_block_size[indI],
												workers[i].cuda_dipole_me,
												workers[i].corr_vector,workers[i].g_threeJ,
												workers[i].g_half_ls[indF][ideg]);
					//cudaDeviceSynchronize();
					//CheckCudaError("K-block");
				}
			}

			double* tmp_half_ls = new double[cuda_dimenmax];

			for(int i = 0; i < workers.size(); i++){
				cudaSetDevice(i);
				cudaDeviceSynchronize();
				cudaMemcpy(tmp_half_ls,workers[i].g_half_ls[indF][ideg],sizeof(double)*cuda_dimenmax,cudaMemcpyDeviceToHost);
				for(int j = 0; j < cuda_dimenmax; j++)
					half_ls[j]+=tmp_half_ls[j];
				
			}		
			for(int i = 0; i < workers.size(); i++)
				cudaMemcpy(workers[i].g_half_ls[indF][ideg],half_ls,sizeof(double)*cuda_dimenmax,cudaMemcpyHostToDevice);
			
			cudaDeviceSynchronize();
			CheckCudaError("Half_ls");
				//cudaMemcpy(half_ls,workers[0].g_half_ls[indF][ideg],sizeof(double)*cuda_dimenmax,cudaMemcpyDeviceToHost);
				//cudaDeviceSynchronize();

}
extern "C" void compute_gpu_half_linestrength(int* indI_,int* indF_,int* ideg_,int* n_,double* vector,double* half_ls){

	



}

/*


extern "C" void compute_gpu_half_linestrength_rotsym_(int* indI_,int* indF_,int* ideg_,double* vector,double* half_ls){

			int indI = *indI_ -1;
			int indF = *indF_ -1;
			int ideg = *ideg_ -1;
			//memcpy(
			cudaMemcpy(corr_vector,vector,sizeof(double)*cuda_dimenmax,cudaMemcpyHostToDevice);
			cudaDeviceSynchronize();
			CheckCudaError("Memcpy");
			int jF = J[indF];
			int jI = J[indI];
			int dJ=jF-jI+1;
			int dimenF = dimen[indF];
			
			if(wigner[indI][dJ]==0) return;
			
			int half_grid_size = (int)ceil((float)((float)dimenF/(float)DIPOLE_BLOCK_SIZE));

			//device_compute_1st_half_ls_flipped_dipole_shared_rotsym_nontrove(
			//const int dimenF,const int dimenI,const int nlevelI,const int nlevelF,const int* kI_,
			//const int* kF_,const int* ktauI_,const int* ktauF_, const int* icorrI_,const int* icorrF_,
			//const double* __restrict__ dipole_me,const double* __restrict__ wigner_,
			//const double* vector,double*  half_ls);

			device_compute_1st_half_ls_flipped_dipole_shared_rotsym_nontrove<<<half_grid_size,DIPOLE_BLOCK_SIZE>>>(dimen[indF],
															dimen[indI],
															nlevelI[indI][dJ],
															nlevelF[indI][dJ],
															K[indI],
															K[indF],
															tau[indI],
															tau[indF],	
															icorr[indI],
															icorr[indF],
															cuda_dipole_me,
															wigner[indI][dJ],
															corr_vector,	
															g_half_ls[indF][ideg]);
			cudaDeviceSynchronize();
			CheckCudaError("Half-Ls");
			cudaMemcpy(half_ls,g_half_ls[indF][ideg],sizeof(double)*cuda_dimenmax,cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();

}
*/

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
	int irootF=blockIdx.x*blockDim.x + threadIdx.x + startF_idx;
	bool valid;

	tauF = 0;
	final_half_ls = 0.0;
	if(irootF<startF_idx+dimenF){

		icontrF =icorrF_[irootF]-1;
		tauF  =  tauF_[irootF] & 1;
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
					final_half_ls+=ls*double(sigmaI)*dipole_me[icontrF + s_icontrI[i]*ncontrF +dipole_idx*dip_stride_2];
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



__global__ void device_compute_1st_half_ls_flipped_dipole_shared_rotsym_nontrove(
const int dimenF,const int dimenI,const int nlevelI,const int nlevelF,const int* kI_,
const int* kF_,const int* ktauI_,const int* ktauF_, const int* icorrI_,const int* icorrF_,
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
	bool valid = false;
	if(irootF<dimenF){
		valid = true;
		icontrF =icorrF_[irootF]-1;
		irlevelF  =  kF_[irootF]-1;
		irdegF   = ktauF_[irootF]-1;
	}



	for(int b_irootI=0; b_irootI < dimenI; b_irootI+=DIPOLE_BLOCK_SIZE){
		irootI = b_irootI+t_id;
		s_icontrI[t_id] = -1;		
		if(irootI < dimenI){
			s_icontrI[t_id]=icorrI_[irootI]-1;
			s_irlevelI[t_id]=kI_[irootI]-1;
			s_irdegI[t_id]=ktauI_[irootI]-1;
			s_vector[t_id]=vector[irootI];
		}
		__syncthreads();
		for(int i = 0; i < DIPOLE_BLOCK_SIZE; i++){
			
			if(s_icontrI[t_id]==-1)continue;
			if(valid){
				icontrI = s_icontrI[i];
				irlevelI= s_irlevelI[i];
				irdegI = s_irdegI[i];
				dip_1 = dipole_me[icontrF + icontrI*dip_stride_1];
				dip_2 = dipole_me[icontrF + icontrI*dip_stride_1 + dip_stride_2];
				dip_3 = dipole_me[icontrF + icontrI*dip_stride_1 + dip_stride_2 + dip_stride_2];
				f_1=wigner_[irlevelI*3 + irlevelF*nlevelI + irdegI*nlevelI*nlevelF + irdegF*g_maxdegen]*dip_1;
				f_2=wigner_[1+irlevelI*3 + irlevelF*nlevelI + irdegI*nlevelI*nlevelF + irdegF*g_maxdegen]*dip_2;
				f_3=wigner_[2+irlevelI*3 + irlevelF*nlevelI + irdegI*nlevelI*nlevelF + irdegF*g_maxdegen]*dip_3;

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


//This relates to TROVEs method of expanding vectors to the primitive basis set
__global__ void device_expand_vectors(const int dimenI,const int igammaI,const int maxcoeff,const int idegI,const int sdeg,const int Ntot,const int* ijterms_, const int* icontr_,const int* N_,const double* repres_, const double* vecI_,double* vec_)
{
    for (int irootI = blockIdx.x * blockDim.x + threadIdx.x; 
         irootI < dimenI; 
         irootI += blockDim.x * gridDim.x) {
	
		int irow,ib,iterm,nelem,isrootI,Ntot,sdeg;
		double dtemp0 = 0.0;
		irow = icontr_[irootI]-1;
		ib = icontr_[irootI + dimenI]-1;
	
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



int main(){}
