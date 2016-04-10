#include "TroveClasses/TroveDipole.h"
#include "TroveClasses/TroveInput.h"
#include "TroveClasses/TroveStates.h"
#include "TroveClasses/TroveBasisSet.h"
#include "BaseClasses/EigenVector.h"
#include "BaseClasses/GpuManager.h"


int main(int argc, char** argv){

	MPI_Init(&argc, &argv);
	TroveDipole test_dipole;
	
	test_dipole.InitDipole(6000000000l);

	TroveInput troveinput;
	troveinput.ReadInput("i_0.1.0.0000.6000.0_A.inp");

	TroveStates troveStates(troveinput);

	troveStates.ReadStates();

	printf("Largest vector is of size: %d\n",troveStates.GetNSizeMax());


	
	GpuManager gpu(0,1);

	TroveBasisSet basisSetj0(0,troveinput.GetNSym(),troveinput.GetSymmetryDegen());
	TroveBasisSet basisSetj1(1,troveinput.GetNSym(),troveinput.GetSymmetryDegen());

	


	basisSetj0.Initialize();
	basisSetj1.Initialize();

	std::vector<int> TryNTot = basisSetj1.GetNTotal();
	for(int i = 0; i < TryNTot.size(); i++){

		printf("%i\n",TryNTot[i]);

	}

	gpu.InitializeAndTransferConstants(troveinput.GetMaxJ(),troveinput.GetNSym(),troveinput.GetSymMaxDegen());
	
	gpu.AllocateVectors(troveinput.GetNJ(),troveStates.GetNSizeMax(),BasisSet::GetDimenMax());

	gpu.TransferBasisSet(&basisSetj0);
	gpu.TransferBasisSet(&basisSetj1);

	gpu.TransferInflation(basisSetj0.GetContr(), basisSetj0.GetIJTerms(),basisSetj0.GetDimensions(),basisSetj0.GetMaxSymCoeffs(),
			basisSetj0.GetMatSize(),basisSetj0.GetRepres(),basisSetj0.GetRepresN(),basisSetj0.GetNTotal(),troveinput.GetSymmetryDegen());

	gpu.TransferInflation(basisSetj1.GetContr(), basisSetj1.GetIJTerms(),basisSetj1.GetDimensions(),basisSetj1.GetMaxSymCoeffs(),
			basisSetj1.GetMatSize(),basisSetj1.GetRepres(),basisSetj1.GetRepresN(),basisSetj1.GetNTotal(),troveinput.GetSymmetryDegen());


	gpu.TransferDipole(&test_dipole,0);


	test_dipole.RemoveDipole();

	EigenVector eigen(troveinput);

	eigen.CacheEigenvectors(&troveStates);


	printf("DimenMax: %i\n",BasisSet::GetDimenMax());
	MPI_Finalize();
	

}
