#pragma once
#include "../BaseClasses/Input.h"
#include "../common/defines.h"
#include "../common/BaseProcess.h"
#include <vector>
class EigenStates{
public:
    int ndeg;
    int igamma;
    int irec;
    int iroot;
    int ilevel;
    double energy;
    int jind;
    int jval;
    int krot;
    int taurot;
    int nsize;
    static bool sort(const EigenStates & lhs, const EigenStates & rhs)
   {
  	return lhs.energy < rhs.energy;
    };
	int rec_len;
    
};


class States : public BaseProcess{
protected:
	int Nmodes;
	std::vector<EigenStates> eigenvalues;
	int Neigenlevels;
	int Neigenroots;
	int sym_maxdegen;
	int sym_nrepres;
	std::vector<int> jVals;
	int minJ;
	int maxJ;
	std::vector<int> sym_gamma;
	std::vector<bool> isym_do;
	std::vector<int> sym_degen;
	std::vector<double> freq_window;
	std::vector<int> igamma_pair;
	std::vector<int> isym_pairs;
	std::vector<double> erange;
	std::vector<double> erange_lower;
	std::vector<double> erange_upper;
	std::vector<double> gns;
	int nsize_max;
	double ZPE;
	bool reduced;
	int symmetry_type;
	bool filter_state(double energy,int igamma);
	bool filter_lower(double energy,int igamma);
	bool filter_upper(double energy,int igamma);
public:
	States(Input & input);
	virtual void ReadStates()=0;
	bool FilterIntensity(int I, int F);
	bool FilterAnyTransitionsFromJ(int I, int J);
	bool FilterLowerState(int I) { return filter_lower(eigenvalues.at(I).energy, eigenvalues.at(I).igamma);}
	bool FilterUpperState(int I) { return filter_upper(eigenvalues.at(I).energy, eigenvalues.at(I).igamma);}
	bool DegeneracyFilter(int gammaI,int gammaF,int idegI,int idegF);
	const char* branch(int jF,int jI);
	int GetNumberStates(){return eigenvalues.size();};
	int GetRecord(int state){return eigenvalues.at(state).irec;};
	int GetRecordLength(int state){return eigenvalues.at(state).rec_len;};
	int GetRoot(int state){return eigenvalues.at(state).iroot;};
	double GetEnergy(int state){return eigenvalues.at(state).energy;};
	int GetJ(int state){return eigenvalues.at(state).jval;};
	int GetJIndex(int state){return eigenvalues.at(state).jind;};
	int GetGamma(int state){return eigenvalues.at(state).igamma;};
	int GetLevel(int state){return eigenvalues.at(state).ilevel;};
	int GetNdeg(int state){return eigenvalues.at(state).ndeg;};
	void GetTransitionDetails(int & num_initial, int & num_trans,int & max_trans);

	int GetNSizeMax(){return nsize_max;};


};
