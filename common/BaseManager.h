#include <cmath>
#include <cstdlib>

#pragma once

class BaseManager{
private:
		size_t available_memory;
		size_t total_memory;
protected:
		void InitializeMemory(size_t bytes);
		void TrackMemory(size_t bytes);
		void FreeMemory(size_t bytes);
		size_t GetAvailableMemory();
		size_t GetAllocatedMemory();
};

