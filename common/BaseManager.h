#include <cmath>
#include <cstdlib>
#include <cstdio>
#pragma once

class BaseManager{
private:
		size_t available_memory;
		size_t total_memory;

		static size_t available_global_memory;
		static size_t total_global_memory;

protected:
		void InitializeMemory(size_t bytes);
		void TrackMemory(size_t bytes);
		void FreeMemory(size_t bytes);
		size_t GetAvailableMemory();
		size_t GetAllocatedMemory();

public:
	static void InitializeGlobalMemory(size_t bytes);
	

	static void TrackGlobalMemory(size_t bytes);
	static void FreeGlobalMemory(size_t bytes);
	static size_t GetAvailableGlobalMemory();
	static size_t GetAllocatedGlobalMemory();


};

