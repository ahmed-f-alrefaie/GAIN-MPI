#include "GpuManager.h"
#include <cmath>

extern "C" double c_three_j(int * j1,int * j2,int * j3,int * k1,int * k2,int * k3);


//Ported trove's three_j subroutine [moltype.f90]
double three_j(int j1,int j2,int j3,int k1,int k2,int k3)
{

	/*
	int newmin,newmax,_new,iphase;
	double a,b,c,al,be,ga,delta,clebsh,minus;
        double term,term1,term2,term3,summ,dnew,term4,term5,term6;
        
      	a = j1;
      	b = j2;
      	c = j3;
      	al= k1;
      	be= k2;
      	ga= k3;

      	double three_j=0.0;
      	
      	if(c > a+b) return three_j;
        if(c < abs(a-b)) return three_j;
        if(a < 0.0 || b < 0.0 || c < 0.0) return three_j;
        if(a < abs(al) || b < abs(be) || c < abs(ga)) return three_j;
        if(-1.0*ga != al+be) return three_j;
      	
      	delta=sqrt(fakt(a+b-c)*fakt(a+c-b)*fakt(b+c-a)/fakt(a+b+c+1.0));     
      	term1=fakt(a+al)*fakt(a-al);
      	term2=fakt(b-be)*fakt(b+be);
      	term3=fakt(c+ga)*fakt(c-ga);
      	term=sqrt((2.0*c+1.0)*term1*term2*term3);
      	
      	newmin=NINT(max(max((a+be-c),(b-c-al)),0.0));
        newmax=NINT(min(min((a-al),(b+be)),(a+b-c))) ;       
        
        summ=0;
        
        for(_new=newmin; _new <=newmax; _new++)
        {
        	dnew=double(_new);
        	term4=fakt(a-al-dnew)*fakt(c-b+al+dnew);
        	term5=fakt(b+be-dnew)*fakt(c-a-be+dnew);
        	term6=fakt(dnew)*fakt(a+b-c-dnew);
        	summ+=pow(-1.0,_new)/(term4*term5*term6);
        }
        
         clebsh=delta*term*summ/sqrt(10.0);
         
         iphase=NINT(a-b-ga);
      	 minus = -1.0;
      	if (iphase%2==0) minus = 1.0;
      	three_j=minus*clebsh/sqrt(2.0*c+1.0);
      	
      	return three_j;
	*/
	return c_three_j(&j1,&j2,&j3,&k1,&k2,&k3);
        
};






GpuManager::GpuManager(int pgpu_id,int nprocs,bool rotsym) : BaseManager(), gpu_id(pgpu_id), Nprocs(nprocs), rotsym_do(rotsym){




	Log("GPUManager Init with ID: %i and Nprocs: %i\n",gpu_id,Nprocs);
	//hls_stream_id = 0;

}


void GpuManager::CheckCudaError(const char* tag){

  cudaSetDevice(gpu_id);
  // check for error
  cudaError_t error = cudaGetLastError();
  if(error != cudaSuccess)
  {
    // print the CUDA error message and exit
    LogErrorAndAbort("[%s] CUDA error: %s\n", tag,cudaGetErrorString(error));
    cudaDeviceReset();
    

  }
};

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
	InitializeMemory(size_t(double(devProp.totalGlobalMem)*0.99));
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


	Log("Transfering ThreeJ symbols");
	//Transfer threeJ symbols
	double* tmp_three_J  = new double[(jmax+1)*(jmax+1)*3*3];

	for(int i =0; i < (jmax+1)*(jmax+1)*3*3; i++)
		tmp_three_J[i]=0.0;

        for(int jI=0; jI <= jmax; jI++)
		for(int jF=std::max(jI-1,0); jF <= std::min(jI+1,jmax); jF++)
			for(int kI=0; kI <= jI; kI++)
				for(int kF=std::max(kI-1,0); kF <= std::min(kI+1,jF); kF++)
				{

					double three = three_j(jI, 1, jF, kI, kF - kI, -kF);
					//printf("(%i,%i,%i,%i) = %14.3E\n",jI,jF,kI,kF,three);
					tmp_three_J[jI + kI*(jmax+1) +(jF-jI + 1)*(jmax+1)*(jmax+1) +  (kF-kI + 1)*(jmax+1)*(jmax+1)*3] = three;
					Log("three-j[%i,%i,%i,%i] = %12.6f\n",jI,jF,kI,kF,tmp_three_J[jI + kI*(jmax+1) +(jF-jI + 1)*(jmax+1)*(jmax+1) +  (kF-kI + 1)*(jmax+1)*(jmax+1)*3]);
				}	
	
	AllocateGpuMemory((void**)&threejsymbols,sizeof(double)*size_t((jmax+1)*(jmax+1)*3*3));
	TransferToGpu(threejsymbols, tmp_three_J,sizeof(double)*size_t((jmax+1)*(jmax+1)*3*3));

	//Transfer symmetry and threeJ constants

	copy_symmetry_constants(SymNrepres,MaxDegen);
	copy_jmax_constant(Jmax);

	Log("Initializing cuBlas..");

	for(int i =0; i < Nprocs; i++){
		handle.push_back(NULL);
		stat = cublasCreate(&handle.back());

		if (stat != CUBLAS_STATUS_SUCCESS) {
			LogErrorAndAbort("CUBLAS initialization failed\n");
		
		}

		cublasSetPointerMode(handle.back(),CUBLAS_POINTER_MODE_DEVICE);
	}
	Log("success!!\n");

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

	AllocateGpuMemory((void**)&tmp_bset.normal_K,sizeof(int)*size_t(tmp_bset.dimensions));
	TransferToGpu(tmp_bset.normal_K, basisSet->GetK(),sizeof(int)*size_t(tmp_bset.dimensions));
	//KBlocks
	if(basisSet->GetKBlock() != NULL){
		AllocateGpuMemory((void**)&tmp_bset.Kblock_size,sizeof(int)*size_t(tmp_bset.J+1));
		TransferToGpu(tmp_bset.Kblock_size, basisSet->GetKBlock(),sizeof(int)*size_t(tmp_bset.J+1));	
	}
	//Vibrational index
	AllocateGpuMemory((void**)&tmp_bset.vib_index,sizeof(int)*size_t(tmp_bset.dimensions));
	TransferToGpu(tmp_bset.vib_index, basisSet->GetVibIndex(),sizeof(int)*size_t(tmp_bset.dimensions));	
	if(basisSet->GetKBlock() != NULL){
	//We need the starting K information here
		tmp_bset.KStart.assign(basisSet->GetKStart(),basisSet->GetKStart()+tmp_bset.J +1);
		tmp_bset.host_KBlock.assign(basisSet->GetKBlock(),basisSet->GetKBlock()+tmp_bset.J +1);
	}

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
		if(rotsym_do)unsorted_half_ls_vectors.push_back(std::vector<double*>());
		for(int deg = 0; deg < MaxDegen; deg++){
			prim_half_ls_vectors.back().push_back(NULL);
			half_ls_vectors.back().push_back(NULL);
			if(rotsym_do)unsorted_half_ls_vectors.back().push_back(NULL);
			AllocateGpuMemory((void**)&prim_half_ls_vectors.back().back(),sizeof(double)*size_t(DimenMax));
			AllocateGpuMemory((void**)&half_ls_vectors.back().back(),sizeof(double)*size_t(DimenMax));
			if(rotsym_do)AllocateGpuMemory((void**)&unsorted_half_ls_vectors.back().back(),sizeof(double)*size_t(DimenMax));
		}
	}

	

	/////////////////final vectors and linestrength
	for(int i = 0; i <Nprocs; i++){
		vectorF.push_back(NULL);
		AllocateGpuMemory((void**)&vectorF.back(),sizeof(double)*size_t(NsizeMax));
		if(rotsym_do){
			unsorted_vectorF.push_back(NULL);
			AllocateGpuMemory((void**)&unsorted_vectorF.back(),sizeof(double)*size_t(DimenMax));
		}
		prim_vectorF.push_back(NULL);
		AllocateGpuMemory((void**)&prim_vectorF.back(),sizeof(double)*size_t(DimenMax));
		linestrength.push_back(NULL);
		AllocateGpuMemory((void**)&linestrength.back(),sizeof(double)*size_t(MaxDegen*MaxDegen));
	}

	//Lets allocate the streams and events here
	for(int i =0; i < nJ; i++){
		half_ls_stream.push_back(std::vector<std::vector<cudaStream_t> >());
		
		correlated_event.push_back(std::vector<cudaEvent_t>());

		completed_half_ls_event.push_back(std::vector<std::vector<cudaEvent_t> >());

		transfer_half_ls_stream.push_back(std::vector<cudaStream_t>());

		correlate_half_ls_stream.push_back(std::vector<cudaStream_t>());

		stream_id.push_back(std::vector<int>());

		for(int deg = 0; deg < MaxDegen; deg++){

			half_ls_stream.back().push_back(std::vector<cudaStream_t>());
			stream_id.back().push_back(0);

			correlated_event.back().push_back(NULL);
			cudaEventCreateWithFlags (&correlated_event.back().back(),cudaEventDisableTiming);

			completed_half_ls_event.back().push_back(std::vector<cudaEvent_t>());
			

			transfer_half_ls_stream.back().push_back(NULL);
			cudaStreamCreate(&transfer_half_ls_stream.back().back());

			correlate_half_ls_stream.back().push_back(NULL);
			cudaStreamCreate(&transfer_half_ls_stream.back().back());

			for(int i = 0; i < MAX_STREAMS; i++){

				half_ls_stream.back().back().push_back(NULL);

				cudaStreamCreate(&half_ls_stream.back().back().back());

				completed_half_ls_event.back().back().push_back(NULL);
				cudaEventCreateWithFlags (&completed_half_ls_event.back().back().back(),cudaEventDisableTiming);
			}
		}
	}
	

	//DotProduct
	for(int i =0; i < Nprocs; i++){
		dot_product_omp_stream.push_back(NULL);
		cudaStreamCreate(&dot_product_omp_stream.back());
	}


	//Dipole streams and events
	cudaEventCreateWithFlags (&transfer_dipole_event,cudaEventDisableTiming);
	cudaStreamCreate(&transfer_dipole_stream);


	Log("Done!!!");
}


void GpuManager::TransferInflation(int* icontr_, int* ijterm,int dimen,int maxsymcoeffs,int matsize,std::vector<double*> repres,std::vector<int*> N,std::vector<int> Ntot,std::vector<int> sym_degen){
	cudaSetDevice(gpu_id);
	GpuInflation tmp_inflate;

	Log("Transferring Inflation\n");
	tmp_inflate.MaxSymCoeffs = maxsymcoeffs;
	
	//Transfer ijtems and icontr
	AllocateGpuMemory((void**)&tmp_inflate.ijTerms,sizeof(int)*tmp_inflate.MaxSymCoeffs*SymNrepres);
	TransferToGpu(tmp_inflate.ijTerms,ijterm,sizeof(int)*tmp_inflate.MaxSymCoeffs*SymNrepres);

	AllocateGpuMemory((void**)&tmp_inflate.contr,sizeof(int)*size_t(dimen)*2l);
	TransferToGpu(tmp_inflate.contr,icontr_,sizeof(int)*size_t(dimen)*2l);


	tmp_inflate.sDeg = sym_degen;

	for(int i = 0; i < SymNrepres; i++){

		tmp_inflate.Ntotal.push_back(Ntot[i]);	
		tmp_inflate.sDeg.push_back(sym_degen[i]);
		tmp_inflate.N.push_back(NULL);

		AllocateGpuMemory((void**)&tmp_inflate.N.back(),sizeof(int)*tmp_inflate.MaxSymCoeffs);	
		TransferToGpu(tmp_inflate.N.back(),N[i],sizeof(int)*tmp_inflate.MaxSymCoeffs);

		

		tmp_inflate.repres.push_back(NULL);
		AllocateGpuMemory((void**)&tmp_inflate.repres.back(),sizeof(double)*Ntot[i]*sym_degen[i]*matsize);	
		TransferToGpu(tmp_inflate.repres.back(),repres[i],sizeof(double)*Ntot[i]*sym_degen[i]*matsize);
	}
	

	
	inflationData.push_back(tmp_inflate);
	Log("done!\n");
	cudaDeviceSynchronize();
	CheckCudaError("Transfer Inflation");


}
/*
void GpuManager::TransferWigner(int dimen,int* root,double* peigenvect){
	cudaSetDevice(gpu_id);
	if(root==NULL || !rotsym_do)
		return;
	
	old_roots.push_back(NULL);
	Log("Sending the old roots and eigenvectros to gpu memory\n");
	AllocateGpuMemory((void**)&old_roots.back(),sizeof(int)*dimen);
	TransferToGpu(old_roots.back(),root,sizeof(int)*dimen);

	eigenvect.push_back(NULL);
	AllocateGpuMemory((void**)&eigenvect.back(),sizeof(double)*dimen);
	TransferToGpu(eigenvect.back(),peigenvect,sizeof(double)*dimen);	
}
*/

void GpuManager::TransferWigner(std::vector<Wigner> p_wigner){
	if(!rotsym_do)
		return;
	Log("Transfering wigner");
	gpu_Wigner.push_back(std::vector<Wigner>());
	for(int i = 0; i < p_wigner.size(); i++){
		Wigner tmp_wigner;
		
		if(p_wigner[i].rot == NULL){
			tmp_wigner.rot = NULL;
			tmp_wigner.nlevelsI = 0;
			tmp_wigner.nlevelsF = 0;
			
		}else{
			tmp_wigner.nlevelsI = p_wigner[i].nlevelsI;
			tmp_wigner.nlevelsF = p_wigner[i].nlevelsF;

			AllocateGpuMemory((void**)&tmp_wigner.rot,sizeof(double)*tmp_wigner.nlevelsI*tmp_wigner.nlevelsF*MaxDegen*MaxDegen*3);
			TransferToGpu(tmp_wigner.rot,p_wigner[i].rot,sizeof(double)*tmp_wigner.nlevelsI*tmp_wigner.nlevelsF*MaxDegen*MaxDegen*3);


		}		
		


		gpu_Wigner.back().push_back(tmp_wigner);





	}



}

void GpuManager::SwitchDipoleBlock(int block){
	cudaSetDevice(gpu_id);
		for(int indF = 0; indF < completed_half_ls_event.size(); indF++)
			for(int idegI = 0; idegI < completed_half_ls_event[indF].size(); idegI++)
				for(int i = 0; i < MAX_STREAMS; i++)
					cudaStreamWaitEvent(transfer_dipole_stream,completed_half_ls_event[indF][idegI][i],0);

			

	
	cudaMemcpyAsync(gpu_dipole,dipole_me->GetDipolePiece(block),dipole_me->GetDipoleSizeBytes(block),cudaMemcpyHostToDevice,transfer_dipole_stream);
	cudaEventRecord(transfer_dipole_event,transfer_dipole_stream);

	dipole_block = block;

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

	//Alloc to the GPU the biggest block size
	
	AllocateGpuMemory((void**)&gpu_dipole,dipole_me->GetBiggestBlockSizeBytes());		
	//Transfer tipole
	TransferToGpu(gpu_dipole,dipole_me->GetDipolePiece(dipole_block),dipole_me->GetDipoleSizeBytes(dipole_block));


	CheckCudaError("Transfer Dipole");
	Log("Done\n");
	
}



void GpuManager::UpdateHalfLinestrength(double* half_ls,int jInd,int ideg){
	cudaSetDevice(gpu_id);
	TransferToGpu(half_ls_vectors[jInd][ideg],half_ls,sizeof(double)*DimenMax);


}	

void GpuManager::UpdateEigenVector(){
	cudaSetDevice(gpu_id);
	TransferToGpu(vectorI,host_vectorI,sizeof(double)*size_t(NsizeMax));
	cudaDeviceSynchronize();
}
void GpuManager::UpdateEigenVector(int proc_id){
	cudaSetDevice(gpu_id);
	cudaMemcpyAsync(vectorF[proc_id],host_vectorF[proc_id],sizeof(double)*size_t(NsizeMax),cudaMemcpyHostToDevice,dot_product_omp_stream[proc_id]);
}

void GpuManager::ExecuteRotSymHalfLs(int indI,int indF,int idegI,int igammaI)
{

	cudaSetDevice(gpu_id);

	int dJ = basisSets[indF].J - basisSets[indI].J + 1;

	if( dJ > 2)
		LogErrorAndAbort("Something went wrong in the dJ department\n");
	
	if(gpu_Wigner[indI][dJ].rot == NULL){
		LogErrorAndAbort("Wigner has a problem\n");
	}
	cudaStreamWaitEvent(half_ls_stream[indF][idegI][0],transfer_dipole_event,0);
	//Execute
	compute_gpu_half_linestrength_rotsym_(
			basisSets[indF].dimensions,
			basisSets[indI].dimensions,

			gpu_Wigner[indI][dJ].nlevelsI,
			gpu_Wigner[indI][dJ].nlevelsF,

			MaxDegen,
			
			basisSets[indI].normal_K,
			basisSets[indF].normal_K,

			basisSets[indI].KTau,
			basisSets[indF].KTau, 

			basisSets[indI].vib_index,
			basisSets[indF].vib_index,
			dipole_me->GetDipoleStart(dipole_block),
			dipole_me->GetDipoleEnd(dipole_block),
			dipole_me->GetDipoleNcontr(dipole_block),
			gpu_dipole,
			
			gpu_Wigner[indI][dJ].rot,
			
			prim_half_ls_vectors[indF][idegI],

			half_ls_vectors[indF][idegI],

			half_ls_stream[indF][idegI][0]);
	cudaEventRecord(completed_half_ls_event[indF][idegI][0],half_ls_stream[indF][idegI][0]);
	cudaStreamWaitEvent(transfer_half_ls_stream[indF][idegI],completed_half_ls_event[indF][idegI][0],0);


	/*cudaDeviceSynchronize();
	double* tmp_half_ls=new double[DimenMax];
	Log("### INDI=%i INDF=%i IGAMMAI=%i#####\n",indI,indF,igammaI); 
	TransferToHost((void*)tmp_half_ls,(void*)half_ls_vectors[indF][idegI],sizeof(double)*basisSets[indF].dimensions);
	cudaDeviceSynchronize();
	for(int i = 0; i < DimenMax; i++)
		Log(" %14.3E\n",tmp_half_ls[i]);
	
	MPI_Abort(MPI_COMM_WORLD,0);
*/
}


void GpuManager::ExecuteHalfLs(int indI,int indF,int idegI,int igammaI){

	TransformHalfLsVector(indI,indF,idegI,igammaI);
	if(rotsym_do)
		ExecuteRotSymHalfLs(indI,indF,idegI,igammaI);
	else
		ExecuteKBlockHalfLs(indI,indF,idegI,igammaI);
	


}


void GpuManager::TransformHalfLsVector(int indI,int indF,int idegI,int igammaI){

	
	//transform to primitive
	transform_vector_primitive(basisSets[indI].dimensions,igammaI,inflationData[indI].MaxSymCoeffs,idegI,inflationData[indI].sDeg[igammaI],inflationData[indI].Ntotal[igammaI],inflationData[indI].ijTerms, inflationData[indI].contr,inflationData[indI].N[igammaI],inflationData[indI].repres[igammaI], vectorI,prim_half_ls_vectors[indF][idegI],correlate_half_ls_stream[indF][idegI]);

	

	cudaEventRecord(correlated_event[indF][idegI],correlate_half_ls_stream[indF][idegI]);


	//Synchronize the streams until the correlation is complete
	for(int i = 0; i < MAX_STREAMS; i++)
		 cudaStreamWaitEvent(half_ls_stream[indF][idegI][i],correlated_event[indF][idegI],0);
	

}



void GpuManager::ExecuteKBlockHalfLs(int indI,int indF,int idegI,int igammaI){
	cudaSetDevice(gpu_id);
	for(int i = 0; i < MAX_STREAMS; i++){
		cudaStreamWaitEvent(half_ls_stream[indF][idegI][i],transfer_dipole_event,0);
	}

	for(int kF=0; kF < basisSets[indF].J+1; kF++){
			compute_gpu_half_linestrength_(
							basisSets[indF].host_KBlock[kF],
							basisSets[indI].dimensions,

							basisSets[indI].J,
							basisSets[indF].J,
	
							kF,
					
							basisSets[indI].KTau,
							basisSets[indF].KTau, 

							basisSets[indI].vib_index,
							basisSets[indF].vib_index,

							basisSets[indF].KStart[kF],
							basisSets[indI].KStart[std::max(kF-1,0)],

							dipole_me->GetDipoleStart(dipole_block),
							dipole_me->GetDipoleEnd(dipole_block),
							dipole_me->GetDipoleNcontr(dipole_block),
							basisSets[indF].host_KBlock[kF],
							basisSets[indI].Kblock_size,

							gpu_dipole,
							prim_half_ls_vectors[indF][idegI],
							threejsymbols,
							half_ls_vectors[indF][idegI],
							half_ls_stream[indF][idegI][GetHStreamId(indF,idegI)]);
		
	}

	for(int i = 0; i < MAX_STREAMS; i++){
		//Record event
		cudaEventRecord(completed_half_ls_event[indF][idegI][i],half_ls_stream[indF][idegI][i]);
		//Push to the parent sream
		cudaStreamWaitEvent(transfer_half_ls_stream[indF][idegI],completed_half_ls_event[indF][idegI][i],0);

	}
	//Now we can Async to the

	/*double* tmp_prim_ls=new double[DimenMax];
	double* tmp_half_ls=new double[DimenMax];
	Log("### INDI=%i INDF=%i IGAMMAI=%i#####\n",indI,indF,igammaI); 
	TransferToHost((void*)tmp_prim_ls,(void*)prim_half_ls_vectors[indF][idegI],sizeof(double)*basisSets[indI].dimensions);
	TransferToHost((void*)tmp_half_ls,(void*)half_ls_vectors[indF][idegI],sizeof(double)*basisSets[indF].dimensions);
	cudaDeviceSynchronize();
	for(int i = 0; i < DimenMax; i++)
		Log("%14.3E -> %14.3E\n",tmp_prim_ls[i],tmp_half_ls[i]);
	*/	
	//MPI_Abort(MPI_COMM_WORLD,0);


}
void GpuManager::ExecuteDotProduct(int indF,int idegI,int idegF,int igammaF,int proc){
	cudaSetDevice(gpu_id);
	//Transform to primitive

		transform_vector_primitive(basisSets[indF].dimensions,igammaF,inflationData[indF].MaxSymCoeffs,idegF,inflationData[indF].sDeg[igammaF],inflationData[indF].Ntotal[igammaF],inflationData[indF].ijTerms, inflationData[indF].contr,inflationData[indF].N[igammaF],inflationData[indF].repres[igammaF], vectorF[proc],prim_vectorF[proc],dot_product_omp_stream[proc]);

	cublasSetStream(handle[proc],dot_product_omp_stream[proc]);

	cublasDdot (handle[proc], basisSets[indF].dimensions,prim_vectorF[proc], 1, half_ls_vectors[indF][idegI], 1, linestrength[proc] + idegI + idegF*MaxDegen );
	//Copy result to host
	


}
