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

