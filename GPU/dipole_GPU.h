
#ifndef DIPOLE_GPU_H_
#define DIPOLE_GPU_H_

#pragma once
#include "cublas_v2.h"
#include <cuda_runtime_api.h>
#include <cuda.h>

void CheckCudaError(const char* tag);


void copy_symmetry_constants(int sym_nrepres,int maxdeg);
void copy_jmax_constant(int jmax);
void copy_dipole_constant(int dip_stride);


//Transforms the symmetry-adapted basis to the primited form
void transform_vector_primitive(const int dimenI,const int igammaI,const int maxcoeff,const int idegI,const int sdeg,const int Ntot,const int* ijterms_, const int* icontr_,const int* N_,const double* repres_, const double* vecI_,double* vec_,cudaStream_t stream);


//Computes a single K block half_linestrength with a certain dipole piece
void compute_gpu_half_linestrength_(const int dimenF,const int dimenI,const int jI,const int jF,
const int kF,const int* tauI_,const int* tauF_, const int* icorrI_,const int* icorrF_,const int startF_idx,int startI_idx,const int startFblock,const int endFblock,const int ncontrF,const int kFBlocksize,const int* kblock_size_,
const double* dipole_me,
const double* vector,const double*  threej,double*  half_ls,cudaStream_t stream);


#endif
