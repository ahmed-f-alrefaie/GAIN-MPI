#include <cstdio>
#include <string>
#include <sys/time.h>
#include <ctime>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <iostream>
#include <vector>
#pragma once

typedef long int int64;
typedef unsigned long int uint64;

void destroy_arr_valid(void** ptr);
inline unsigned int Fortran2D_to_1D(int i, int j,int isize){ return i+ j*isize;};
inline unsigned int Fortran3D_to_1D(int i, int j,int k,int isize,int jsize){ return i+ j*isize + k*isize*jsize ;};
size_t GetFilenameSize(std::string name);
bool fexists(const char *filename);
void CreateRandomVector(double** vector,size_t count);
// trim from start
std::string &ltrim(std::string &s);

// trim from end
std::string &rtrim(std::string &s);

// trim from both ends
std::string &trim(std::string &s);

std::vector<std::string> split(std::string const &input);

void ReadFortranRecord(FILE* IO, void* data_ptr);
double readd();
int readi();
char* readc();
void assertdouble(double & d1, double & d2, double tol);
double GenRandNumber(double LO, double HI);
int64 GetTimeMs64();
void CreateZeroVector(double** vector,size_t count);
void wrapvalue(unsigned long int & var, unsigned long int min,unsigned long int max);
