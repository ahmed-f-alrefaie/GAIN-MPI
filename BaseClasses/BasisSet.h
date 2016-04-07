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
    int 	dimen;
    int* 	ktau;
    int* 	k;
    int*	vib_index;
    int* 	ktau;
    int* 	kstart;
    static TO_BasisSet* j0Basis;

    
public:
   TO_BasisSet(int pid,int max_pros,int J);
   virtual void Initialize() = 0;

   const int GetJval(){return jval;};
   const int GetMaxDegSize(){return max_deg_size;};
   const int GetDimensions(){return Maxcontracts;};
   const int * vib_index(){return iroot_correlat_j0;};
   const int * GetK(){return k;};
};
