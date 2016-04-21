#include "MultiGpuManager.h"
#include <climits>
#pragma once

void MultiGpuManager::PinVectorMemory(double* vector,int n){
	cudaHostRegister(vector,size_t(n)*sizeof(double),cudaHostRegisterPortable);
}

void MultiGpuManager::UnpinVectorMemory(double* vector,int n){
	cudaHostUnregister(vector);
}

int MultiGpuManager::GetFreeDevice(){
	int device_id=-1;

	int devCount;
	cudaGetDeviceCount(&devCount);
	for(device_num; device_num< devCount; device_num++){
		cudaSetDevice(device_num );
		if(cudaFree(0)==cudaSuccess){
			cudaThreadExit();
			return device_num++;
		}	
	}

	return -1;


}

void MultiGpuManager::InitializeAndTransferConstants(int jmax,int sym_repres,int pmax_degen){
	DegenMax = pmax_degen;
	for(int i = 0; i < m_gpus.size(); i++)
		m_gpus[i]->InitializeAndTransferConstants(jmax,sym_repres,pmax_degen);	

}



MultiGpuManager::MultiGpuManager(std::vector<int> jvals,States* states,int nprocs,bool rotsym) : BaseProcess(){

	total_gpus_assigned = 0;
	int allowed_gpus = 0;
	device_num = 0;
	const char* num_gpu = getenv("NUM_GPUS");

	if(num_gpu!= NULL){
		allowed_gpus = atoi(num_gpu);		

	}
	
	if(allowed_gpus <= 0){
		
	}
	
	int gpu_id = 0;
	while(gpu_id!=-1){
		gpu_id = GetFreeDevice();
		if(gpu_id!=-1){
			//LogErrorAndAbort("Error could not get device!\n");
			m_gpus.push_back(new GpuManager(gpu_id,nprocs,rotsym));
			total_gpus_assigned++;
		}
		if(total_gpus_assigned == allowed_gpus && (allowed_gpus > 0)){
			break;
		}
	}



	VerboseLog("We have found and will utilize %d GPUs\n",total_gpus_assigned);

	m_states = states;
	m_jvals = jvals;
}

void MultiGpuManager::TransferBasisSet(BasisSet* basisSet){
	for(int i = 0; i < m_gpus.size(); i++)
		m_gpus[i]->TransferBasisSet(basisSet);


}
void MultiGpuManager::TransferInflation(int* icontr_, int* ijterm,int dimen,int maxsymcoeffs,int matsize,std::vector<double*> repres,std::vector<int*> N,std::vector<int> Ntot,std::vector<int> sym_degen){

	for(int i = 0; i < m_gpus.size(); i++)
		m_gpus[i]->TransferInflation(icontr_, ijterm,dimen,maxsymcoeffs,matsize,repres,N,Ntot,sym_degen);

}
void MultiGpuManager::TransferDipole(Dipole* dipole_)
{


	m_dipole = dipole_;
	
	size_t smallest_memory = ULLONG_MAX;

	for(int i = 0; i < total_gpus_assigned; i++){
		smallest_memory = std::min(m_gpus[i]->GetAvailableMemory(),smallest_memory);

	}
	
	Log("The smallest amount of memory we have is %12.6f GB\n",float(smallest_memory)/float(GIGABYTE_TO_BYTE)); 
	
	Log("Initializing dipole.......\n");
	m_dipole->InitDipole(smallest_memory);

	int cur_gpu=0;
	//Lets distrubute it evenly
	for(int i = 0; i < m_dipole->GetNumBlocks(); i++){
		if(cur_gpu >= total_gpus_assigned) cur_gpu=0;		
		m_dipole_dist.push_back(cur_gpu++);
		Log("Dipole block %d is being given to the %d(st/th) GPU we are managing\n",i,cur_gpu-1);

	}
	cur_gpu=0;
	//Now lets give each GPU heir repective blocks
	for(int i = 0; i < m_dipole->GetNumBlocks(); i++){
		if(cur_gpu >= total_gpus_assigned) break;
		m_gpus[cur_gpu++]->TransferDipole(m_dipole,i);
	}
	//If its not blocked or we could put everything into memory then clean up the dipole
	if(!m_dipole->IsBlocked() || m_dipole_dist.size() <= total_gpus_assigned){
		m_dipole->RemoveDipole();
	}else{
		//Lets pin the memory for further transfers
		for(int i = 0; i < m_dipole->GetNumBlocks(); i++){
			PinVectorMemory(m_dipole->GetDipolePiece(i),m_dipole->GetDipoleSizeBytes(i));
		}

	}
	Log("Done\n");





}

void MultiGpuManager::AllocateVectors(int nJ,int nsizemax,int dimenmax){
	for(int i = 0; i < m_gpus.size(); i++)
		m_gpus[i]->AllocateVectors(nJ, nsizemax,dimenmax);

	NsizeMax = nsizemax;
	DimenMax = dimenmax;
	this->nJ = nJ;
	cudaMallocHost((void**)&vectorI,sizeof(double)*size_t(nsizemax));

	for(int i = 0; i < nJ; i++){
		tmp_half_linestrength.push_back(std::vector<double*>());
		half_linestrength.push_back(std::vector<double*>());
		for(int idegI = 0; idegI < DegenMax; idegI++){
			half_linestrength.back().push_back(NULL);
			cudaMallocHost((void**)&half_linestrength.back().back(),size_t(DimenMax)*sizeof(double));
			tmp_half_linestrength.back().push_back(NULL);
			cudaMallocHost((void**)&tmp_half_linestrength.back().back(),size_t(DimenMax)*sizeof(double));
			//PinVectorMemory(half_linestrength.back().back(),size_t(DimenMax)*sizeof(double));
			//PinVectorMemory(tmp_half_linestrength.back().back(),size_t(DimenMax)*sizeof(double));
			BaseManager::TrackGlobalMemory(sizeof(double)*DimenMax);
			BaseManager::TrackGlobalMemory(sizeof(double)*DimenMax);
			
		}
	}	


}
	
void MultiGpuManager::UpdateHalfLinestrength(double* half_ls,int jInd,int ideg)
{
	for(int i = 0; i < m_gpus.size(); i++)
		m_gpus[i]->UpdateHalfLinestrength(half_ls,jInd,ideg);


}

double* MultiGpuManager::GetHalfLineStrength(int indF,int idegI)
{
	return half_linestrength.at(indF).at(idegI);
}
	
double* MultiGpuManager::GetInitialVector(){
	return vectorI;
}
double* MultiGpuManager::GetFinalVector(int proc_id){
	int gpu_id = proc_id % total_gpus_assigned;
	return m_gpus[gpu_id]->GetFinalVector(proc_id);

}
double* MultiGpuManager::GetLinestrength(int proc_id){
	int gpu_id = proc_id % total_gpus_assigned;
	return m_gpus[gpu_id]->GetLinestrength(proc_id);

}


void MultiGpuManager::UpdateEigenVector(){
	for(int i = 0; i < m_gpus.size(); i++){
		memcpy((void*)m_gpus[i]->GetInitialVector(),vectorI,sizeof(double)*size_t(NsizeMax));
		m_gpus[i]->UpdateEigenVector();
	}
}
void MultiGpuManager::UpdateEigenVector(int proc_id){
	int gpu_id = proc_id % total_gpus_assigned;
	m_gpus[gpu_id]->UpdateEigenVector(proc_id);

}

void MultiGpuManager::TransferWigner(std::vector<Wigner> p_wigner){
	for(int i = 0; i < m_gpus.size(); i++)
		m_gpus[i]->TransferWigner(p_wigner);
}

void MultiGpuManager::ExecuteHalfLs(int iLevelI,int indI,int ndegI,int igammaI,int igammaF)
{
	//transform the vector anyway
	for(int i = 0; i < m_gpus.size(); i++)
		for(int indF = 0; indF < nJ; indF++)
			for(int idegI = 0; idegI < ndegI; idegI++)
				m_gpus[i]->TransformHalfLsVector(indI,indF,idegI,igammaI);

	//Lets loop through the Dipole blocks
	for(int i = 0; i < m_dipole->GetNumBlocks(); i++){
		int selected_gpu = m_dipole_dist[i];
		//Lets wait for the GPU to finish first
		m_gpus[selected_gpu]->WaitForDevice();
		//If it has the right dipole block then we push it to the device
		if(m_gpus[selected_gpu]->GetCurrentBlock() != i) m_gpus[i]->SwitchDipoleBlock(i);
		//Otherwise do the calculations
		for(int indF = 0; indF < nJ; indF++){
			if(!m_states->FilterAnyTransitionsFromJ(iLevelI,m_jvals[indF]))
				continue;
			
			for(int idegI = 0; idegI < ndegI; idegI++){
				//Execute the half linestrength
				m_gpus[selected_gpu]->ExecuteHalfLs(indI,indF, idegI,igammaI);
			}
		}
	}
	for(int indF = 0; indF < nJ; indF++){
		//Get the result from zero
		for(int idegI = 0; idegI < ndegI; idegI++){
			m_gpus[0]->GetHalfLineStrengthResult(half_linestrength[indF][idegI],indF,idegI);
		}
	}

	//Get the half linestrength result and apply it to the real half linestrength, we skip this if only one gpu has the result
	for(int i = 1; i < total_gpus_assigned; i++){
		for(int indF = 0; indF < nJ; indF++){
		//Get the result from zero
			for(int idegI = 0; idegI < ndegI; idegI++){
				m_gpus[i]->GetHalfLineStrengthResult(tmp_half_linestrength[indF][idegI],indF,idegI);
				#pragma omp parallel for
				for(int n=0; n< DimenMax; n++)
					half_linestrength[indF][idegI][n]+=tmp_half_linestrength[indF][idegI][n];
			}
		}		

	}



}
void MultiGpuManager::ExecuteDotProduct(int indF,int idegI,int idegF,int igammaF,int proc)
{
	int gpu_id = proc % total_gpus_assigned;
	m_gpus[gpu_id]->ExecuteDotProduct(indF,idegI,idegF,igammaF,proc);	

}
void MultiGpuManager::WaitForLineStrengthResult(int proc_id)
{
	int gpu_id = proc_id % total_gpus_assigned;
	m_gpus[gpu_id]->WaitForLineStrengthResult(proc_id);

}




