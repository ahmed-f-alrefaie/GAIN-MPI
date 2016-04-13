#pragma once

#include "../common/BaseProcess.h"
#include "../BaseClasses/States.h"
#include "../common/defines.h"



class Output  : public BaseProcess {


private:
	char filename[1024];
	FILE* output;
	int maxDeg;
	States* states;	
	bool to_file;
	double threshold;
	double Acoef_s_1;
	double Q;
	double temperature;
	double ZPE;
	bool reduced;
	double beta;

public:
	Output(States* pstates,double ptemperature,double ppartition,double thresh,int pmaxDeg,bool red,const char* pfilename=NULL);
	
	void Initialize();

	void OutputLinestrength(int ilevelI,int ilevelF,double* linestrength);

	void Close();

	~Output();
		





};
