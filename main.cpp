#include "TroveClasses/TroveDipole.h"
#include "TroveClasses/TroveInput.h"
#include "TroveClasses/TroveStates.h"
#include "TroveClasses/TroveBasisSet.h"
#include "BaseClasses/EigenVector.h"
#include "BaseClasses/GpuManager.h"
#include "BaseClasses/Output.h"
#include "common/Timer.h"
#include <omp.h>
#include <vector>
#include <cstring>
#include <string>

int main(int argc, char** argv){

	MPI_Init(&argc, &argv);

	int rank;
	int nProcs;
	bool quit_prog = false;
	bool do_file = false;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&nProcs);	
	char input_filename[1024];
	char output_filename[1024];

	int omp_threads = omp_get_max_threads();
	int nJ;

	printf("OMP_NUM_THREADS=%i\n",omp_threads);

	/////////////////////Simple argument handling///////////////

	//Lets read the arguments
	if(rank ==0){
		if(argc<2){
			printf("Usage is [-f] <input file> [<output filename>]\n");
			quit_prog = true;

		}else{
			if(argc==2){
				strcpy(input_filename,argv[1]);
			}else if(argc==4){
				if(std::string(argv[1]) == "-f"){
					strcpy(input_filename,argv[2]);	
					strcpy(output_filename,argv[3]);
					do_file = true;	
				}else{
					printf("Usage is [-f] <input file> [<output filename>]\n");
					quit_prog = true;
				}
					
			}else{
				printf("Usage is [-f] <input file> [<output filename>]\n");
				quit_prog = true;
			}
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
	Output* m_output;



	//Read inputs
	m_input = new TroveInput();
	m_input->ReadInput(input_filename);

	nJ = m_input->GetNJ();

	//ReadStates
	
	m_states = new TroveStates((*m_input));
	m_states->ReadStates();


	//Initialize output
	if(do_file)
		m_output = new Output(m_states,m_input->GetTemperature(),m_input->GetPartition(),m_input->GetThreshold(),m_input->GetSymMaxDegen(),output_filename);
	else
		m_output = new Output(m_states,m_input->GetTemperature(),m_input->GetPartition(),m_input->GetThreshold(),m_input->GetSymMaxDegen());

	m_output->Initialize();
	std::vector<int> m_jvals = m_input->GetJvals();
	GpuManager* m_gpu;
	for(int i=0; i < nProcs; i++){
	//Initialize our gpu
		if(rank==i){
			m_gpu = new GpuManager(-1,omp_threads);
			m_gpu->InitializeAndTransferConstants(m_input->GetMaxJ(),m_input->GetNSym(),m_input->GetSymMaxDegen());
		}
		MPI_Barrier( MPI_COMM_WORLD);
	}
	
	

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

	int DimenMax = BasisSet::GetDimenMax();


	//Alocate needed vectors
	m_gpu->AllocateVectors(m_input->GetNJ(),m_states->GetNSizeMax(),BasisSet::GetDimenMax());


	//Handle dipole
	m_dipole = new TroveDipole();
	m_dipole->InitDipole(6000000000l);

	m_gpu->TransferDipole(m_dipole,0);

	m_dipole->RemoveDipole();



	int num_initial;
	int num_trans;
	int max_trans;
	

	m_states->GetTransitionDetails(num_initial, num_trans,max_trans);

	
	int nLevels = m_states->GetNumberStates();

	double* vector_I = m_gpu->GetInitialVector();
	

	std::vector< std::vector<double*> > half_linestrength;

	for(int i = 0; i < nJ; i++){

		half_linestrength.push_back(std::vector<double*>());

		for(int idegI = 0; idegI < m_input->GetSymMaxDegen(); idegI++){
			half_linestrength.back().push_back(NULL);
			half_linestrength.back().back() = new double[DimenMax];
			BaseManager::TrackGlobalMemory(sizeof(double)*BasisSet::GetDimenMax());
			
		}
	}

	eigen = new EigenVector((*m_input));

	eigen->CacheEigenvectors(m_states);

	double ZPE = m_input->GetZPE();


	MPI_Barrier( MPI_COMM_WORLD);
	if(rank ==0){
		printf("ZPE: %12.6f \n ",ZPE);

		printf("-------------------------------------Begin Intensity Calculation----------------------------\n");
	}
	MPI_Barrier( MPI_COMM_WORLD);
	int transitions = 0;
	
	for(int iLevelI = 0; iLevelI < nLevels; iLevelI++){
		
		//All MPI processes should do this
		if(!m_states->FilterLowerState(iLevelI))
			continue;

		int jI= m_states->GetJ(iLevelI);
		int gammaI = m_states->GetGamma(iLevelI);
		int nSizeI = m_states->GetRecordLength(iLevelI);
		int ndegI = m_states->GetNdeg(iLevelI); 
		int indI =  m_states->GetJIndex(iLevelI);
		double energyI = m_states->GetEnergy(iLevelI);
		int indexI = m_states->GetLevel(iLevelI);
		int expected_process = eigen->ReadVector(vector_I,iLevelI,nSizeI);
		
		m_gpu->UpdateEigenVector();


		if(rank==0) {printf("Lower state energy: %12.6f cm-1\n",energyI-ZPE); fflush(0);}

		Timer::getInstance().StartTimer("Half linestrength");
		for(int indF = 0; indF < nJ; indF++){

			if(!m_states->FilterAnyTransitionsFromJ(iLevelI,m_jvals[indF]))
				continue;

			for(int idegI = 0; idegI < ndegI; idegI++){			
				if(expected_process == rank){
					//Do halflinestrength
					m_gpu->ExecuteHalfLs(half_linestrength[indF][idegI],indI,indF,idegI,gammaI);
				
				}
				
				

			

				//Broadcast
				MPI_Bcast(half_linestrength[indF][idegI], DimenMax, MPI_DOUBLE, expected_process, MPI_COMM_WORLD);
				//Update our gpu
				
				if(expected_process != rank) m_gpu->UpdateHalfLinestrength(half_linestrength[indF][idegI],indF,idegI);


				
			

			}
			//MPI_Abort(MPI_COMM_WORLD,0);

		}

	
		Timer::getInstance().EndTimer("Half linestrength");

		Timer::getInstance().StartTimer("Intensity Loop");
		
		#pragma omp parallel for reduction(+:transitions)
		for(int iLevelF=rank; iLevelF < nLevels; iLevelF+=nProcs){

			if(!m_states->FilterIntensity(iLevelI,iLevelF))
				continue;
	
			
			int thread_id = omp_get_thread_num();

			int jF= m_states->GetJ(iLevelF);
			int gammaF = m_states->GetGamma(iLevelF);
			int nSizeF = m_states->GetRecordLength(iLevelF);
			int ndegF = m_states->GetNdeg(iLevelF); 
			int indF =  m_states->GetJIndex(iLevelF);
			double energyF = m_states->GetEnergy(iLevelF);
			int indexF = m_states->GetLevel(iLevelF);


			double* vector_F = m_gpu->GetFinalVector(thread_id);
			double* linestrength = m_gpu->GetLinestrength(thread_id);

			if(rank != eigen->ReadVector(vector_F,iLevelF,nSizeF)){
				printf("Error!!!!\n");
				exit(0);
			}
		

			//Do work
			//Output result
			m_gpu->UpdateEigenVector(thread_id);



			for(int idegF = 0; idegF < ndegF; idegF++)
				for(int idegI = 0; idegI < ndegI; idegI++)
					m_gpu->ExecuteDotProduct(indF,idegI,idegF,gammaF,thread_id);

			m_gpu->WaitForLineStrengthResult(thread_id);

			double linestr,ls; 
			ls = 0.0;

			m_output->OutputLinestrength(iLevelI,iLevelF,linestrength);

			transitions++;

			

		}
		MPI_Barrier( MPI_COMM_WORLD);
		Timer::getInstance().EndTimer("Intensity Loop");
		fflush(0);
		int g_transitions =0;
		MPI_Reduce(&transitions,&g_transitions,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
		
		if(rank==0){
			float total_time_hours = 1.0/3600.0;
			
			float lines_per_second = float(g_transitions)/Timer::getInstance().GetTotalTimeInSeconds("Intensity Loop");
			
			float predicted_time = 	(float(num_trans)/lines_per_second)*total_time_hours;
			float percentage_of_transitions = (float)transitions*100.0/(float)g_transitions;
			printf("### %d / %d transitions, L/s: %14.3E Pred Total time %8.4f hours [%12.6f]\n",g_transitions,num_trans,lines_per_second,predicted_time,percentage_of_transitions);					
					
			
		}


		MPI_Barrier( MPI_COMM_WORLD);
		

	}
	m_output->Close();
	MPI_Barrier( MPI_COMM_WORLD);
	if(rank==0){
		Timer::getInstance().PrintTimerInfo();
	}

	MPI_Finalize();
	

}
