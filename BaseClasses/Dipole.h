#include "../common/defines.h"
#include "../common/BaseProcess.h"
#include "../common/BaseManager.h"
#include <vector>
#include <cmath>
#pragma once


struct DipolePiece{
	double* dipole_me;
	int startF;
	int endF;
	int ncontrF;
	size_t size;
};


//Reads the dipole_ flipped
class Dipole : public BaseProcess {
protected:
	std::vector<DipolePiece> dipole_me;
	size_t dipole_size;
	int MaxContracts;
	int num_blocks;
	size_t max_size;
	bool block_form;
public:		
	Dipole(size_t avail_mem);

	virtual void InitDipole()=0;
	//const double* GetDipole(){ return dipole_me;}
	size_t GetDipoleSize(){ return dipole_size;}
	size_t GetDipoleSizeBytes(){return dipole_size*sizeof(double);}
	int GetNumBlocks(){return num_blocks;};
	bool IsBlocked(){return block_form;}
	const double* GetDipolePiece(int block){ return dipole_me.at(block).dipole_me;};
	int GetDipoleStart(int block){ return dipole_me.at(block).startF;};
	int GetDipoleEnd(int block){ return dipole_me.at(block).endF;};
	int GetDipoleNcontr(int block){ return dipole_me.at(block).ncontrF;};
	size_t GetDipoleSize(int block){ return dipole_me.at(block).size;};
	size_t GetMaxContracts(){ return MaxContracts;};
	virtual void CleanUp();
	~Dipole();

};
