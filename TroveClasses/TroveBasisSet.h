#pragma once
#include "../common/defines.h"
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <fstream>
#include "../common/Util.h"
#include "../BaseClasses/BasisSet.h"


//#define DEBUG

struct TO_PTrotquantaT
{
     int	 j;  
     int	 k;
     int	tau;
};
struct TO_PTintcoeffsT 
{      
      
      char*	type;
      int*	icoeffs;  //iCoeffs indexes - arbitrary information
      int	size1,size2;
      
};

struct TO_PTrepresT
{
	double* repres;
	int* N;
};



class TroveBasisSet : public BasisSet
{

private:

    int 	Maxsymcoeffs;
    int 	Nclasses;
    int ncases;
    int nlambdas;
    int ncontr;
    int nclasses;
    int icontr;
    int iroot;

    int* 	icontr2icase;
    int*	icase2icontr;
    int*	contractive_space;
    TO_PTintcoeffsT*  index_deg;
    TO_PTrotquantaT* rot_index;
    int*	icontr_correlat_j0;
    int* 	nsize;
    int nsize_max;
    TO_PTrepresT* irr;
    int* Ntotal;
    int mat_size;
    int* ijterms;
    //Required for correlation
    static TroveBasisSet* j0BasisSet;

    std::vector<int> sym_degen;
    int sym_nrepres;
    
    
    //for later when i might change the basis set
    virtual void Correlate();
    void ComputeijTerms();
    void ReadBasisSet();
public:
   TroveBasisSet(int J,int sym_n, std::vector<int> psym_degen);

   void Initialize();
   
   int GetMaxSymCoeffs(){
	return Maxsymcoeffs;
   }

   int GetMatSize(){return mat_size;};
   int* GetContr(){return icontr2icase; };

   int* GetIJTerms(){return ijterms; };

   std::vector<int> GetSymDegen(){return sym_degen;};
   std::vector<int> GetNTotal(){ 
	std::vector<int> NtotTmp;
	NtotTmp.assign(Ntotal,Ntotal+sym_nrepres);
	return NtotTmp;
  };


  std::vector<double*> GetRepres(){
	std::vector<double*> RepresTmp;
	for(int i =0; i < sym_nrepres; i++)
		RepresTmp.push_back(irr[i].repres);

	return RepresTmp;
	};
  std::vector<int*> GetRepresN(){
	std::vector<int*> RepresTmp;
	for(int i =0; i < sym_nrepres; i++)
		RepresTmp.push_back(irr[i].N);

	return RepresTmp;
	};
};



