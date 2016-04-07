#include "TroveDipole.h"
#include "../common/Util.h"
#include <cstdio>
#include <cstring>

TroveDipole::TroveDipole(size_t mem) : Dipole(mem){}

//Reads the dipole_ flipped
void TroveDipole::InitDipole(){
	char buff[20];
	int imu,imu_t;
	int ncontr_t;
   	size_t matsize,rootsize,rootsize_t,rootsize2;
	FILE* extF = fopen("j0_extfield.chk","r");

	if(extF == NULL)
	{
		Log("[read_dipole] Error: checkpoint file of name %s not found","j0_extfield.chk");
		//fprintf(stderr,"[read_dipole] Error: checkpoint file not found");
		exit(0);
	}
	
	ReadFortranRecord(extF, buff);
	
	if(memcmp(buff,"Start external field",20) != 0)
	{
		Log("[read_dipole] Error: checkpoint file of name %s has bogus header %s","j0_extfield.chk");
		//fprintf(stderr,"[read_dipole] bogus header");
		exit(0);
	}
	
	ReadFortranRecord(extF, &ncontr_t);
	
	/////////////////////////
	
	
   	rootsize2 = ncontr_t*ncontr_t;
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
			Log("[read_dipole] has bogus imu - restore_vib_matrix_elements: %i /= %i",imu_t,(i+1));
			//fprintf(stderr,"[read_dipole] bogus imu");
			exit(0);
		}
		
		ReadFortranRecord(extF, (temp_dipole) + i*rootsize2);
				
			 
	}
		
	ReadFortranRecord(extF, buff);		 
		if(memcmp(buff,"End external field",18) != 0)
	{
		printf("[read_dipole] Error: checkpoint file of name %s has bogus footer %s","j0_extfield.chk",buff);
		fprintf(stderr,"[read_dipole] bogus footer");
		exit(0);
	}
	
	fclose(extF);


	#ifdef DEBUG
	//for(int i = 0; i < ncontr_t; i++)
//		for(int j = 0; j < ncontr_t; j++)//
			//for(int k = 0; k < 3; k++)
		//		printf("dipole[%i,%i,%i] = %16.8e\n",i,j,k,(*dipole_me)[i + j*ncontr_t + k*ncontr_t*ncontr_t]);
	//exit(0);
	#endif
	int n_contr_block = ceil(float(ncontr_t)/float(num_blocks));
	int cur_block_size = 0;
	int startF=0,endF=0,ncontrF=0;
	//Allocate X parts
	//dipole_block.dip_block=new FDipole_block[parts];
	//dipole_block.parts=parts;	
	Log("Flipping dipole...and blocking into %i pieces\n",num_blocks);
	for(int blocks = 0; blocks < num_blocks; blocks++){
		startF=blocks*n_contr_block;
		
		endF = (blocks+1)*n_contr_block;
		endF = std::min(endF,ncontr_t);
		ncontrF = endF-startF;
		dipole_me.push_back(DipolePiece());
		//Allocate the dipole
		dipole_me.back().startF = startF;
		dipole_me.back().endF = endF;
		dipole_me.back().dipole_me = new double[ncontrF*ncontr_t*3];
		dipole_me.back().size = sizeof(double)*size_t(ncontrF)*size_t(ncontr_t)*3l;
		BaseManager::TrackGlobalMemory(dipole_me.back().size);
		Log("Size of block %i is %zu\n",blocks,dipole_me.back().size);
		dipole_me.back().ncontrF = ncontrF;
		Log("startF =%i, block number = %i, n_contr_block = %i\n",startF,blocks,ncontrF);
		for(int i = 0; i < ncontr_t; i++)
			for(int f = startF; f < endF; f++)
				for(int k = 0; k < 3; k++)
					dipole_me.back().dipole_me[f-startF + i*ncontrF + k*ncontr_t*ncontrF] = temp_dipole[i + f*ncontr_t + k*ncontr_t*ncontr_t]; 

	}


	delete[] temp_dipole;
	printf("done!\n");


};



