#include "EigenVector.h"

EigenVector::EigenVector(Input & input) : BaseProcess(){
	
	jVals = input.GetJvals();
	isym_do=input.GetISymDo();
	sym_nrepres = input.GetNSym();


}
void EigenVector::CacheEigenvectors(States* pstates){
	//Lets allocate memory in the heap;
	//Lets get the states
	Log("Caching Eigenvectors...");
	states = pstates;
	
	//Log("allocating heap......done!\n");
	
	Log("Total memory for eigenvectors is %12.6f GB\n",((double)BaseManager::GetAvailableGlobalMemory())/1e9);

	//total_vals = BaseManager::GetAvailableGlobalMemory()/sizeof(double);
	cur_vals = 0;
	
	char filename[1024];

	for(int j = 0; j < jVals.size(); j++){
		//Push back a J file
		eigenvector_files.push_back(std::vector<MPI_File>());
		for(int g = 0; g < sym_nrepres; g++){

			eigenvector_files.back().push_back(NULL);
			if(!isym_do[g])
				continue;
			
			//We open it
			sprintf(filename,"j0eigen_vectors%d_%d.chk",jVals[j],g+1);
			Log("Opening %s\n",filename);
			
			if(MPI_File_open(MPI_COMM_WORLD,filename,MPI_MODE_RDONLY,MPI_INFO_NULL,&eigenvector_files.back().back()) != MPI_SUCCESS){
				LogErrorAndAbort("Could not open %s!!\n",filename);
				
			}

		}


	}


	//vector_heap = new double[total_vals];

	//double * heap_ptr = vector_heap;
	cached_vectors = 0;
	total_vectors=0;
	//total_vectors = states->GetNumberStates()/m_num_processes;
	//count vectors
	for(int i = m_process_id; i < states->GetNumberStates(); i+=m_num_processes)
		total_vectors++;

	//Let us begin
	for(int i = m_process_id; i < states->GetNumberStates(); i+=m_num_processes){

		//if (i % m_num_processes != m_process_id)
		//	continue;
		//Lets read it		
		int jInd = states->GetJIndex(i);
		int gamma = states->GetGamma(i);
		int record = states->GetRecord(i);
		int rec_len = states->GetRecordLength(i);
		//
		//total_vectors++;
		//Can we fit it?
		//if(cur_vals + rec_len > total_vals)
		//	break;

		//Otherwise we can continue
		size_t memory_needed = size_t(rec_len)*sizeof(double);

		if(memory_needed >= BaseManager::GetAvailableGlobalMemory())
			break;


		stored_vectors.push_back(NULL);
		stored_vectors.back() = new double[rec_len];
		BaseManager::TrackGlobalMemory(size_t(rec_len)*sizeof(double));		
		
		MPI_Status status;
		MPI_File_read_at(eigenvector_files[jInd][gamma],size_t(record)*size_t(rec_len)*sizeof(double),stored_vectors.back(),size_t(rec_len),MPI_DOUBLE,&status);
		//Read
		//fseek(eigenvector_files[jInd][gamma],size_t(record)*size_t(rec_len)*sizeof(double),SEEK_SET);
		
		//fread(stored_vectors.back(),sizeof(double),size_t(rec_len),eigenvector_files[jInd][gamma]);

		cached_vectors++;
		
	}
	Log("We have cached %d out of %d vectors here!\n",cached_vectors,total_vectors);
	if(cached_vectors==total_vectors){
		Log("Alright! we've cached all vectors here...Lets close the files. Hopefully all other processes did the same! Lets check...\n");
	}

	int global_vector_count;

	MPI_Reduce(&cached_vectors,&global_vector_count,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	Log("Vector count: %d\n",global_vector_count);
	if(global_vector_count == states->GetNumberStates()){
		Log("Yes!!! We've cached everything! No more IO! Lets close everything\n");
		for(int j = 0; j < jVals.size(); j++){
			for(int g = 0; g < sym_nrepres; g++){
				if(!isym_do[g])
					continue;
				if(eigenvector_files[j][g]==NULL)
					continue;
				MPI_File_close(&eigenvector_files[j][g]);
			}
			
		}
		eigenvector_files.clear();

	}

}

void EigenVector::ReadVectorFromFile(double* array,int nLevel,size_t size){

		int jInd = states->GetJIndex(nLevel);
		int gamma = states->GetGamma(nLevel);
		int record = states->GetRecord(nLevel);
		int rec_len = states->GetRecordLength(nLevel);
		//Otherwise we can continue
		//fseek(eigenvector_files[jInd][gamma],size_t(record)*size_t(rec_len)*sizeof(double),SEEK_SET);
		MPI_Status status;
		MPI_File_read_at(eigenvector_files[jInd][gamma],size_t(record)*size_t(rec_len)*sizeof(double),array,size_t(rec_len),MPI_DOUBLE,&status);
		//fread(array,sizeof(double),size_t(rec_len),eigenvector_files[jInd][gamma]);	

		


}

void EigenVector::ReadVectorFromHeap(double* array,int nLevel,size_t size){

	memcpy(array,stored_vectors.at(nLevel),size*sizeof(double));

}


int EigenVector::ReadVector(double* array,int nLevel,size_t size){
	
	int which_proc = nLevel % m_num_processes;

	if (which_proc  == m_process_id){
		
		int level=nLevel/m_num_processes;

		if(level < stored_vectors.size()){
			
			//Read from  cahce
			ReadVectorFromHeap(array,level,size);
			
		
		}else{
			ReadVectorFromFile(array,nLevel,size);
		
		}
//

	}

	return which_proc;



}



