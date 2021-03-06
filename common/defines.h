#pragma once


#include <mpi.h>

#define GIGABYTE_TO_BYTE 1000000000l
#define MEGABYTE_TO_BYTE 1000000l

#ifdef KEPLER
	#define MAX_STREAMS 32
#else
	#define MAX_STREAMS 16
#endif
//Math constants

//#define DEBUG
#define PI 3.14159265359
#define AVOGNO 6.0221415e+23
#define PLANCK 6.62606896e-27 
#define VELLGT 2.99792458e+10
#define BOLTZ 1.380658e-16
#define ACOEF 3.13618923e-7
#define SQRT2 0.70710678118
#define CMCOEF 1.327209365e-12
