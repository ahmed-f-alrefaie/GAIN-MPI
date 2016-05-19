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
	bool output_linestrengths;
	bool compute_intensities;
	double threshold;
	std::vector<double> gns;
	double Q;
	double cmcoef;
	double temperature;
	double ZPE;
	bool reduced;
	double beta;


public:
	Output(States* pstates,std::vector<double> p_gns,double ptemperature,double pZPE,double ppartition,double thresh,int pmaxDeg,bool red,bool do_file=false,const char* pfilename=NULL,bool full_line=false,bool compute_intens=false);
	
	void Initialize();

	void OutputLinestrength(int ilevelI,int ilevelF,double* linestrength);

	void Close();

	~Output();
		





};
