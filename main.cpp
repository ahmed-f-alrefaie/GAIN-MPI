#include "TroveClasses/TroveDipole.h"

int main(int argc, char** argv){

	MPI_Init(&argc, &argv);
	TroveDipole test_dipole(4000000000l);
	
	test_dipole.InitDipole();

	MPI_Finalize();
	

}
