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
	
	Log("allocating heap......");
	
	total_vals = BaseManager::GetAvailableGlobalMemory()/sizeof(double);
	cur_vals = 0;
	
	char filename[1024];

	for(int j = 0; j < jVals.size(); j++){
		//Push back a J file
		eigenvector_files.push_back(std::vector<FILE*>());
		for(int g = 0; g < sym_nrepres; g++){

			eigenvector_files.back().push_back(NULL);
			if(!isym_do[g])
				continue;
			Log("Opening files for %d,%d\n",j,g);
			//We open it
			sprintf(filename,"j0eigen_vectors%d_%d.chk",jVals[j],g+1);
			eigenvector_files.back().back() = fopen(filename,"rb");
			if(eigenvector_files.back().back() == NULL){
				Log("Could not open %s!!\n",filename);
				exit(0);
			}

		}


	}


	vector_heap = new double[total_vals];

	double * heap_ptr = vector_heap;
	cached_vectors = 0;
	total_vectors=0;
	//total_vectors = states->GetNumberStates()/m_num_processes;

	//Let us begin
	for(int i = 0; i < states->GetNumberStates(); i++){

		if (i % m_num_processes != m_process_id)
			continue;
		//Lets read it		
		int jInd = states->GetJIndex(i);
		int gamma = states->GetGamma(i);
		int record = states->GetRecord(i);
		int rec_len = states->GetRecordLength(i);
		total_vectors++;
		//Can we fit it?
		if(cur_vals + rec_len > total_vals)
			break;

		//Otherwise we can continue

		//Read
		fseek(eigenvector_files[jInd][gamma],record,SEEK_SET);
		
		fread(heap_ptr,sizeof(double),size_t(rec_len),eigenvector_files[jInd][gamma]);
		//Put the information
		vector_heap_information.push_back({heap_ptr,rec_len});

		//Increment the pointer
		heap_ptr+=rec_len;
		//Count the cache
		cached_vectors++;
	}
	Log("We have cached %d out of %d vectors here!\n",cached_vectors,total_vectors);
	if(cached_vectors==total_vectors){
		Log("Alright! we've cached all vectors here...Hopefully all other processes did the same! Lets check...");
	}

	int global_vector_count;

	MPI_Reduce(&cached_vectors,&global_vector_count,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
	Log("Vector count: %d\n",global_vector_count);
	if(global_vector_count == states->GetNumberStates()){
		Log("Yes!!! We've cached everything! No more IO!\n");


	}

}

bool EigenVector::ReadVectorFromFile(double* array,int nLevel,size_t size){

		int jInd = states->GetJIndex(nLevel);
		int gamma = states->GetGamma(nLevel);
		int record = states->GetRecord(nLevel);
		int rec_len = states->GetRecordLength(nLevel);
		//Otherwise we can continue
		fseek(eigenvector_files[jInd][gamma],record,SEEK_SET);
		
		fread(array,sizeof(double),size_t(rec_len),eigenvector_files[jInd][gamma]);	




}



bool EigenVector::ReadVector(double* array,int nLevel,size_t size){
	
	if (nLevel % m_num_processes != m_process_id){
		return false;

	}else{

		int level=nLevel/m_num_processes;

		if(level < cached_vectors){
			//Read from file
			return ReadVectorFromFile(array,nLevel,size);
		
		}else{
			if(size != vector_heap_information[level].vector_size){
				Log("Problem, size discrepency!\n");
				exit(0);
			}
			memcpy(array,vector_heap_information[level].start,size_t(vector_heap_information[level].vector_size)*sizeof(double));
			return true;
		
		}


	}

}



