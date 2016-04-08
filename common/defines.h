#pragma once

#include <mpi.h>


//Kepler definitions
#ifdef KEPLER
	#define CORRELATE_BLOCK_SIZE 512
	#define DIPOLE_BLOCK_SIZE 128
	#define MAX_STREAMS 32
#else
	#define CORRELATE_BLOCK_SIZE 512
	#define DIPOLE_BLOCK_SIZE 128
	#define MAX_STREAMS 16
#endif

#define GIGABYTE_TO_BYTE 1000000000l
#define MEGABYTE_TO_BYTE 1000000l


//Math constants

//#define DEBUG

#define WARP_SIZE 32
#define PI 3.14159265359
#define AVOGNO 6.0221415e+23
#define PLANCK 6.62606896e-27
#define VELLGT 2.99792458e+10
#define BOLTZ 1.380658e-16
#define ACOEF 3.13618923e-7
#define SQRT2 0.70710678118
