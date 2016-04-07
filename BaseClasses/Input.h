#include "../common/defines.h"
#include "../common/BaseProcess.h"
#include "../common/BaseManager.h"
#include "../common/Util.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>

#pragma once
class Input : public BaseProcess {
protected:
	double ZPE;
	std::vector<double> erange;
	std::vector<double> erange_lower;
	std::vector<double> erange_upper;
	std::vector<double>freq_window;
	std::vector<int> jVals;
	std::vector<bool> isym_do;
	std::vector<int> isym_pairs;
	double q_stat;
	double temperature;
	std::vector<double> gns;
	//Fmole_type molec;
	std::string filename;
	std::string output_file;
	double thresh_linestrength;
	bool reduced;
	int nJ;	
	std::string symmetry_group;
	std::ifstream stream;
	std::vector<int> igamma_pair;
	int nmodes;


	

	//Symmetry information
	int maxJ;
	int sym_nrepres;
	std::vector<int> sym_degen;
	int sym_maxdegen;
	std::vector<std::string> c_sym;



	//Size
	size_t memory;

public:
	Input() : BaseProcess(){};
	virtual void ReadInput(const char* filename)=0;
	
	//Getters
	double GetZPE(){return ZPE;};

	size_t GetGlobalMemory(){return memory;};

	std::vector<double> GetErange(){return erange;};
	std::vector<double> GetErangeLower(){return erange_lower;};
	std::vector<double> GetErangeUpper(){return erange_upper;};
	std::vector<double> GetFreqWindow(){return freq_window;};
	std::vector<int> GetJvals(){return jVals;};
	double GetNmodes(){return nmodes;};
	std::vector<int> GetSymmetryDegen(){return sym_degen;};
	int GetNSym(){ return sym_nrepres;};
	int GetSymMaxDegen(){return sym_maxdegen;} ;
	std::vector<bool> GetISymDo(){return isym_do;};
	std::vector<int> GetIGammaPair(){return igamma_pair;};
	std::vector<int> GetISymPairs(){return isym_pairs;};
	//std::vector< std::vector<int> > GetQuantaLower(){ return quanta_lower;}
	//std::vector< std::vector<int> > GetQuantaUpper(){ return quanta_upper;}
	int GetMaxJ(){ return maxJ;}
	int GetNJ(){return nJ;};
	double GetTemperature(){return temperature;};
	double GetPartition(){return q_stat;};
	double GetGNS(int gamma){return gns.at(gamma);};
	std::vector<double> GetGNS(){return gns;};
	bool DoSym(int gamma){return isym_do.at(gamma);}
	double GetThreshold(){return thresh_linestrength;}
	

};
