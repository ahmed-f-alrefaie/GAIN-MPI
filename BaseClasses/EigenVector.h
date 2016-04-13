#include <map>
#include "../common/BaseProcess.h"
#include "../common/BaseManager.h"
#include "../BaseClasses/States.h"
#include "../common/defines.h"
#include <cstdio>
#include "../common/Util.h"
#pragma once


struct vector_info{
 	double* start;
	size_t vector_size;
};

class EigenVector : public BaseProcess {
private:

	int max_vector_count_mem;
	int cached_vectors;
	int total_vectors;
	size_t total_vals;
	size_t cur_vals;

	//std::vector<double> vectors
	//double* vector_heap;
	std::vector<double*> stored_vectors;
	States* states;
	//std::vector<vector_info> vector_heap_information;
	std::vector< std::vector<MPI_File> > eigenvector_files;
		
	std::vector<int> jVals;
	std::vector<bool> isym_do;
	int sym_nrepres;
	void ReadVectorFromFile(double* array,int nLevel,size_t size);
	virtual void ReadVectorFromHeap(double* array,int nLevel,size_t size);
	//
public:
	EigenVector(Input & input);
	virtual void CacheEigenvectors(States* pstates);
	int ReadVector(double* array,int nLevel,size_t size);
	void Close();


};
