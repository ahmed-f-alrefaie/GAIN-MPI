#include "TroveClasses/TroveDipole.h"
#include "TroveClasses/TroveInput.h"
#include "TroveClasses/TroveStates.h"
#include "TroveClasses/TroveBasisSet.h"
#include "BaseClasses/EigenVector.h"
#include "BaseClasses/GpuManager.h"
#include <omp.h>
#include <vector>
#include <cstring>
int main(int argc, char** argv){

	MPI_Init(&argc, &argv);

	int rank;
	bool quit_prog = false;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	char input_filename[1024];

	int omp_threads = omp_get_max_threads();

	printf("OMP_NUM_THREADS=%i\n",omp_threads);

	/////////////////////Simple argument handling///////////////


	//Lets read the arguments
	if(rank ==0){
		if(argc<2){
			printf("Usage is <input file>\n");
			quit_prog = true;

		}else{
			strcpy(input_filename,argv[1]);
		}




	}


	MPI_Bcast(&quit_prog, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
	MPI_Bcast(input_filename, 1024, MPI_CHAR, 0, MPI_COMM_WORLD);

	if(quit_prog){
		MPI_Finalize();
		return 0;
	}

	Input* m_input;
	States* m_states;
	std::vector<TroveBasisSet*> m_basisSets;
	Dipole* m_dipole;
	EigenVector* eigen;


	


	//Read inputs
	m_input = new TroveInput();
	m_input->ReadInput(input_filename);

	//ReadStates
	
	m_states = new TroveStates((*m_input));
	m_states->ReadStates();

	std::vector<int> m_jvals = m_input->GetJvals();
	
	//Initialize our gpu
	GpuManager* m_gpu = new GpuManager(0,omp_threads);

	
	m_gpu->InitializeAndTransferConstants(m_input->GetMaxJ(),m_input->GetNSym(),m_input->GetSymMaxDegen());

	for(int i = 0; i < m_input->GetNJ(); i++){
		m_basisSets.push_back(new TroveBasisSet(m_jvals[i],m_input->GetNSym(),m_input->GetSymmetryDegen()));
		m_basisSets.back()->Initialize();
		//Transfer to the gpu
		m_gpu->TransferBasisSet(m_basisSets.back());
		//Transfer inflation
		TroveBasisSet* t_b = (TroveBasisSet*)m_basisSets.back();
		m_gpu->TransferInflation(t_b->GetContr(), t_b->GetIJTerms(),t_b->GetDimensions(),t_b->GetMaxSymCoeffs(),
			t_b->GetMatSize(),t_b->GetRepres(),t_b->GetRepresN(),t_b->GetNTotal(),m_input->GetSymmetryDegen());

		
	}
	//Alocate needed vectors
	m_gpu->AllocateVectors(m_input->GetNJ(),m_states->GetNSizeMax(),BasisSet::GetDimenMax());


	//Handle dipole
	m_dipole = new TroveDipole();
	m_dipole->InitDipole(6000000000l);

	m_gpu->TransferDipole(m_dipole,0);

	m_dipole->RemoveDipole();

	eigen = new EigenVector((*m_input));

	eigen->CacheEigenvectors(m_states);

	








/*

	//Read the input file
	TroveInput troveinput;
	troveinput.ReadInput(input_filename);

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

*/
	MPI_Finalize();
	

}
