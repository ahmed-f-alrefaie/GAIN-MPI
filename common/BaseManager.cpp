#include "BaseManager.h"
#include <cstdlib>

size_t BaseManager::available_global_memory;
size_t BaseManager::total_global_memory;

void BaseManager::InitializeMemory(size_t bytes){
	available_memory = bytes;
	total_memory = bytes;
}

void BaseManager::TrackMemory(size_t bytes){
	
	available_memory -= bytes;
	if(available_memory < 0){
		printf("No more memory!!\n");
	}


}
void BaseManager::FreeMemory(size_t bytes){
	available_memory += bytes;
}
size_t BaseManager::GetAvailableMemory(){
	return available_memory;
}
size_t BaseManager::GetAllocatedMemory(){
	return total_memory-available_memory;
}


void BaseManager::InitializeGlobalMemory(size_t bytes){
		BaseManager::available_global_memory = bytes;
		BaseManager::total_global_memory = bytes;
}

void BaseManager::TrackGlobalMemory(size_t bytes){
	
	BaseManager::available_global_memory -= bytes;
		if(BaseManager::available_global_memory < 0){
			printf("No more memory!!\n");
		}



}
void BaseManager::FreeGlobalMemory(size_t bytes){
		BaseManager::available_global_memory += bytes;
}
size_t BaseManager::GetAvailableGlobalMemory(){
		return BaseManager::available_global_memory;
}
size_t BaseManager::GetAllocatedGlobalMemory(){
		return BaseManager::total_global_memory-BaseManager::available_global_memory;
}
