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
	size_t GetBiggestBlock();
	bool block_form;
public:		
	Dipole();

	virtual void InitDipole(size_t avail_mem)=0;
	//const double* GetDipole(){ return dipole_me;}
	size_t GetBiggestBlockSize() { 
		size_t t_size=0;
		for(int i = 0; i < dipole_me.size(); i++){
			t_size = std::max(t_size,dipole_me[i].size);
		}
		return t_size;

	};
	size_t GetBiggestBlockSizeBytes(){return GetBiggestBlockSize()*sizeof(double);};
	size_t GetDipoleSize(){ return dipole_size;}
	size_t GetDipoleSizeBytes(){return dipole_size*sizeof(double);}
	int GetNumBlocks(){return num_blocks;};
	bool IsBlocked(){return (num_blocks > 1);}
	double* GetDipolePiece(int block){ return dipole_me.at(block).dipole_me;};
	int GetDipoleStart(int block){ return dipole_me.at(block).startF;};
	int GetDipoleEnd(int block){ return dipole_me.at(block).endF;};
	int GetDipoleNcontr(int block){ return dipole_me.at(block).ncontrF; };
	size_t GetDipoleSize(int block){ return dipole_me.at(block).size; };
	size_t GetDipoleSizeBytes(int block){ return GetDipoleSize(block)*sizeof(double); };
	size_t GetMaxContracts(){ return MaxContracts;};
	void RemoveDipole();
	~Dipole();

};
