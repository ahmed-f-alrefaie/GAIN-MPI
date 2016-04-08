#pragma once
#include "../common/defines.h"
#include <iostream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <fstream>
#include "../common/Util.h"
#include "../common/BaseProcess.h"

class BasisSet : public BaseProcess
{

protected:
    int 	jval;
    int 	max_deg_size;
    int 	Maxcontracts;
    int* 	ktau;
    int* 	k;
    int*	iroot_correlat_j0;
    int* 	kstart;
    int*        k_block_size; 
    static int dimenMax;
public:
   BasisSet(int J);
   virtual void Initialize() = 0;

   int GetJval(){return jval;};
   int GetMaxDegSize(){return max_deg_size;};
   int GetDimensions(){return Maxcontracts;};
   const int * GetVibIndex(){return iroot_correlat_j0;};
   const int * GetK(){return k;};
   const int * GetKTau(){return ktau;};
   const int * GetKBlock(){return k_block_size;};
   const int * GetKStart(){return kstart;};
   static int GetDimenMax(){return dimenMax;};
};
