#include "TroveClasses/TroveDipole.h"
#include "TroveClasses/TroveInput.h"
#include "TroveClasses/TroveStates.h"
#include "TroveClasses/TroveBasisSet.h"
#include "BaseClasses/EigenVector.h"
int main(int argc, char** argv){

	MPI_Init(&argc, &argv);
	TroveDipole test_dipole(6000000000l);
	
	test_dipole.InitDipole();

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
	MPI_Finalize();
	

}
