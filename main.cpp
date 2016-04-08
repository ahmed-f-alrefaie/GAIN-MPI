#include "TroveClasses/TroveDipole.h"
#include "TroveClasses/TroveInput.h"
#include "TroveClasses/TroveStates.h"
#include "TroveClasses/TroveBasisSet.h"
#include "BaseClasses/EigenVector.h"
int main(int argc, char** argv){

	MPI_Init(&argc, &argv);
	TroveDipole test_dipole;
	
	test_dipole.InitDipole(6000000000l);

	TroveInput troveinput;
	troveinput.ReadInput("i_0.1.0.0000.6000.0_A.inp");

	TroveStates troveStates(troveinput);

	troveStates.ReadStates();

	EigenVector eigen(troveinput);

	eigen.CacheEigenvectors(&troveStates);
	

	TroveBasisSet basisSetj0(0,troveinput.GetNSym(),troveinput.GetSymmetryDegen());
	TroveBasisSet basisSetj1(1,troveinput.GetNSym(),troveinput.GetSymmetryDegen());

	


	basisSetj0.Initialize();
	basisSetj1.Initialize();

	std::vector<int> TryNTot = basisSetj1.GetNtotal();
	for(int i = 0; i < TryNTot.size(); i++){

		printf("%i\n",TryNTot[i]);

	}

	printf("DimenMax: %i\n",BasisSet::GetDimenMax());
	MPI_Finalize();
	

}
