#include <map>
#include "../common/BaseProcess.h"
#include "../common/defines.h"
#include <cstdio>
#include "../common/Util.h"
#pragma once


struct vector_info{
 	double* start;
	size_t vector_size;
}

class EigenVector : public BaseProcess, public BaseManager {
private:

	int max_vector_count_mem;
	//std::vector<double> vectors
	double* vector_heap;
	std::vector<vector_info> vector_heap_information;
	//
public:
	EigenVector(size_t available_memory);
	void CacheEigenvectors(States* states);
	void ReadVector(double* array,int J, int G, int record,size_t size);


};
