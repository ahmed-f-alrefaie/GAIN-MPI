#include "Dipole.h"
#pragma once



Dipole::Dipole(){
	block_form = false;
}


void Dipole::RemoveDipole(){
	Log("Cleaning up dipole for extra memory\n");
	for(int i = 0; i < dipole_me.size(); i++){
		delete[] dipole_me[i].dipole_me;
		BaseManager::FreeGlobalMemory(dipole_me[i].size);
	}

}



Dipole::~Dipole(){}
