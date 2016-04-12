#include "TroveStates.h"
#include <string>
#include<fstream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include "../common/Util.h"

const char* j0eigen_descr_gamma_filebase="j0eigen_descr%d_%d.chk";

TroveStates::TroveStates(Input & pInput) : States(pInput){}


void TroveStates::ReadJGStates(int jVal,int gamma,int jInd){
	
	char filename[1024];
	std::string line;
	//The firs process is to count how many states pass the initial tests
	int nlevels = 0;
	int nroots = 0;
	int nroots_t=0;	
	int dim_basis=0;
	char* line_ptr;
	int igamma,ideg,maxdeg,irec, ilevel,ilarge_coef;
	double energy;
	int maxj = jVal;
	int iroot = 0;
	maxdeg = 0;


	//std::vector<int> ilevel_new;	
	int (*ktau_rot)[2] = new int[1 + (2*maxj)][2];

	for(int k0 = 1; k0 <= jVal; k0++){
		ktau_rot[0][0] = 0;
		ktau_rot[0][1] = jVal % 2;       
	//ktau_rot(0,1) = 0
        //ktau_rot(0,2) = mod(Jval(jind),2)
		iroot = 0;
		/*
			do k0 = 1,Jval(jind)
		*/
		for(int k0= 1; k0 <= jVal; k0++)
		{
			for(int tau0 = 0; tau0 <=1; tau0++){
			   iroot = iroot + 1;
			   ktau_rot[iroot][0] = k0;
			   ktau_rot[iroot][1] = tau0;
			   
			}
		}
	}

	sprintf(filename,j0eigen_descr_gamma_filebase,jVal,(gamma+1));
	std::ifstream descr_file(filename);
	if(!descr_file)
	{
		LogErrorAndAbort("Error! couldn';t open %s!!!!\n",filename);
		
	}
	Log("I'm sorry but I can't read the fingerprints of %s yet so I'm skipping them :( please be careful ;_;\n",filename);
	while(trim(line).compare("Start Quantum numbers and energies")!=0)
	{
		getline(descr_file,line);
		if(descr_file.eof()){
			LogErrorAndAbort("Error! malformed descr file %s!!!!\n",filename);
			
		}
					
	}
	getline(descr_file,line);
			//Get the nroots
	getline(descr_file,line);
	getline(descr_file,line);
	nroots_t = strtol(line.c_str(),&line_ptr,0);
	dim_basis = strtol(line_ptr,&line_ptr,0);
	nsize_max = std::max(nsize_max,dim_basis);
	Log("J=%i nroots = %i, dim_basis = %i\n",jVal,nroots_t,dim_basis);


	while(getline(descr_file,line)){

				//If we hit the signiture then we are done
				if(trim(line).compare("End Quantum numbers and energies")==0 || descr_file.eof())
					break;
	
				int irec = strtol(line.c_str(),&line_ptr,0)-1;//irec
				int igamma = strtol(line_ptr,&line_ptr,0)-1; //igamma
				if((igamma) != gamma){
					Log("Gammas dont match\n");
					exit(0);
				}
				int ilevel = strtol(line_ptr,&line_ptr,0)-1; //ilevel
				int ideg=strtol(line_ptr,&line_ptr,0)-1; //ideg
				//Get the energy
				double energy = strtod(line_ptr,&line_ptr);
				if(!filter_state(energy,gamma))
					continue;

				int k_ind = strtol(line_ptr,&line_ptr,0);

				//Log("Nmodes %d MaxDegen %d\n",Nmodes,sym_maxdegen);
				eigenvalues.push_back(EigenStates());
				eigenvalues.back().irec = irec;
				eigenvalues.back().jind       = jInd;
				eigenvalues.back().jval       = jVal;
				eigenvalues.back().ilevel     = ilevel;
				eigenvalues.back().energy     = energy;
				eigenvalues.back().igamma     = gamma;
				eigenvalues.back().rec_len = dim_basis;

				eigenvalues.back().krot       = ktau_rot[k_ind][0];
				eigenvalues.back().taurot     = ktau_rot[k_ind][1];
				eigenvalues.back().ndeg  = sym_degen.at(gamma);	
					
		
	}

					
	descr_file.close();					

	delete[] ktau_rot;
};


void TroveStates::ReadStates(){

	//Read states
	for(int i=0; i < jVals.size(); i++){
		int J = jVals[i];
		for(int g = 0; g < sym_nrepres; g++){
			if(isym_do[g]){
				ReadJGStates(J,g,i);
			}
			
		}



	}



	if(eigenvalues.size()==0){
		LogErrorAndAbort("Filters are too tight!\n");
		
	}

	Log("Sorting %d states.................",GetNumberStates());
	std::sort(eigenvalues.begin(),eigenvalues.end(),EigenStates::sort);
	Log("..done\n");

	Neigenlevels=eigenvalues.size();

}


