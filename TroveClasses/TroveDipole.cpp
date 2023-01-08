#include "TroveDipole.h"
#include "../common/Util.h"
#include <cstdio>
#include <cstring>

TroveDipole::TroveDipole() : Dipole(){}

//Reads the dipole_ flipped
void TroveDipole::InitDipole(size_t avail_mem){

	max_size = avail_mem;

	char buff[20];
	size_t imu,imu_t;
	size_t ncontr_t;
   	size_t matsize,rootsize,rootsize_t,rootsize2;
	FILE* extF = fopen("j0_extfield.chk","r");

	if(extF == NULL)
	{
		LogErrorAndAbort("[read_dipole] Error: checkpoint file of name %s not found","j0_extfield.chk");
		//fprintf(stderr,"[read_dipole] Error: checkpoint file not found");
		
	}
	
	ReadFortranRecord(extF, buff);
	
	if(memcmp(buff,"Start external field",20) != 0)
	{
		LogErrorAndAbort("[read_dipole] Error: checkpoint file of name %s has bogus header %s","j0_extfield.chk");
		//fprintf(stderr,"[read_dipole] bogus header");
		
	}
	
	ReadFortranRecord(extF, &ncontr_t);
	
	/////////////////////////
	
	MaxContracts = ncontr_t;

   	rootsize2 = size_t(ncontr_t)*size_t(ncontr_t);
	matsize = rootsize2*3;
	dipole_size = matsize;
	size_t total_mat_bytes = matsize*sizeof(double);
	num_blocks = (total_mat_bytes/max_size) + 1;
	rootsize  = ncontr_t*(ncontr_t+1)/num_blocks;
	Log("We need to split the matrix into %d blocks\n",num_blocks);
	double* temp_dipole = new double[matsize];
	//(*dipole_me) = new double[matsize];
	int contr_div = ceil(float(ncontr_t)/float(num_blocks));	
	for(int i = 0 ; i< 3; i++)
	{
		ReadFortranRecord(extF, &imu_t);
		if(imu_t != (i+1))
		{
			LogErrorAndAbort("[read_dipole] has bogus imu - restore_vib_matrix_elements: %i /= %i",imu_t,(i+1));
			//fprintf(stderr,"[read_dipole] bogus imu");
			
		}
		
		size_t read_records = ReadFortranRecord(extF, (temp_dipole) + i*rootsize2);
		Log("Expected %zu bytes, actually read %zu bytes\n",rootsize2*8l,read_records);
				
			 
	}
		
	ReadFortranRecord(extF, buff);		 
		if(memcmp(buff,"End external field",18) != 0)
	{
		LogErrorAndAbort("[read_dipole] Error: checkpoint file of name %s has bogus footer %s","j0_extfield.chk",buff);
	}
	
	fclose(extF);


	#ifdef DEBUG
	//for(int i = 0; i < ncontr_t; i++)
//		for(int j = 0; j < ncontr_t; j++)//
			//for(int k = 0; k < 3; k++)
		//		printf("dipole[%i,%i,%i] = %16.8e\n",i,j,k,(*dipole_me)[i + j*ncontr_t + k*ncontr_t*ncontr_t]);
	//exit(0);
	#endif
	size_t n_contr_block = ceil(float(ncontr_t)/float(num_blocks));
	int cur_block_size = 0;
	size_t startF=0,endF=0,ncontrF=0;

	//Trove creates dipoles in reversed order so we need to flip them
	Log("Flipping dipole...and blocking into %i pieces\n",num_blocks);
	for(size_t blocks = 0; blocks < num_blocks; blocks++){
		startF=blocks*n_contr_block;
		
		endF = (blocks+1)*n_contr_block;
		endF = std::min(endF,size_t(ncontr_t));
		ncontrF = endF-startF;
		dipole_me.push_back(DipolePiece());
		//Allocate the dipole
		dipole_me.back().startF = startF;
		dipole_me.back().endF = endF;
		dipole_me.back().dipole_me = new double[ncontrF*ncontr_t*3];
		dipole_me.back().size = size_t(ncontrF)*size_t(ncontr_t)*3l;
		BaseManager::TrackGlobalMemory(dipole_me.back().size);
		Log("Size of block %i is %zu\n",blocks,dipole_me.back().size);
		dipole_me.back().ncontrF = ncontrF;
		Log("startF =%i, block number = %i, n_contr_block = %i\n",startF,blocks,ncontrF);
		for(size_t i = 0; i < ncontr_t; i++)
			for(size_t f = startF; f < endF; f++)
				for(size_t k = 0; k < 3; k++)
					dipole_me.back().dipole_me[f-startF + i*ncontrF + k*ncontr_t*ncontrF] = temp_dipole[i + f*ncontr_t + k*ncontr_t*ncontr_t]; 

	}


	delete[] temp_dipole;
	Log("done!\n");


};



