#pragma once
#include "../common/defines.h"
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <fstream>
#include "../common/Util.h"
#include "../BaseClasses/BasisSet.h"
#include "../common/Wigner.h"

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

class K_tau_vib{
public:
	int K;
	int tau;
	int vib;
	int iroot;
	int dimen;
	double eigenvec;
	static bool sortK(const K_tau_vib & lhs, const K_tau_vib & rhs)
   	{
  		return lhs.K < rhs.K || ((lhs.K==rhs.K) && lhs.iroot< rhs.iroot);
    	};
	static bool sortRoot(const K_tau_vib & lhs, const K_tau_vib & rhs)
   	{
  		return lhs.iroot < rhs.iroot;
    	};
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

    int* 	icontr2icase;
    int*	icase2icontr;
    int*	contractive_space;
    TO_PTintcoeffsT*  index_deg;
    TO_PTrotquantaT* rot_index;
    int*	icontr_correlat_j0;
    int* 	nsize;
    int*        original_root;
    int nsize_max;
    TO_PTrepresT* irr;
    int* Ntotal;
    int mat_size;
    int* ijterms;
    bool do_rotsym;
    //Required for correlation
    static TroveBasisSet* j0BasisSet;

    std::vector<int> sym_degen;
    int sym_nrepres;
    double* eigenvects;
    int* dimenRoot;
    
    std::vector<Wigner> wigner;


    int* old_roots;
    //for later when i might change the basis set
    virtual void Correlate();
    void ComputeijTerms();
    void ReadBasisSet();
public:
   TroveBasisSet(int J,int sym_n, std::vector<int> psym_degen,bool rotsym=false);

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

  int* GetOldRoots(){if(do_rotsym)return old_roots;else return NULL;};
  double* GetRotEigenvects(){if(do_rotsym)return eigenvects;else return NULL;};

  std::vector<Wigner> GetWigner(){return wigner;};

};



