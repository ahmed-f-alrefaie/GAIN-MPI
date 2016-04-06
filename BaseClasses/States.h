class EigenStates{
public:
    int ndeg;
    int igamma;
    std::vector<int> irec;
    std::vector<int> iroot;
    int ilevel;
    double energy;
    int jind;
    int jval;
    int krot;
    int taurot;
    EigenStates(int maxDegen){
	irec.resize(maxDegen);
	iroot.resize(maxDegen);

    };
    static bool sort(const EigenStates & lhs, const EigenStates & rhs)
   {
  	return lhs.energy < rhs.energy;
    };
	int rec_len;
    
};


class States : public BaseProcess{
private:
	int Nmodes;
	std::vector<Eigen> eigenvalues;
	int Neigenlevels;
	int Neigenroots;
	int sym_maxdegen;
	int Nclasses;
	std::vector<int> jVals;
	int minJ;
	int maxJ;
	std::vector<int> sym_gamma;
	std::vector<bool> isym_do;
	std::vector<int> sym_degen;
	std::vector<double> freq_window;
	std::vector<int> igamma_pair;
	std::vector<int> isym_pairs;
	std::vector<std::string> c_sym;
	std::vector<double> erange;
	std::vector<double> erange_lower;
	std::vector<double> erange_upper;
	std::vector< std::vector<int> >quanta_lower;
	std::vector< std::vector<int> > quanta_upper;	
	std::vector<double> gns;
	int nsize_max;
	double ZPE;
	bool filter(double energy,int igamma);
public:
	States(Input & input);
	virtual void ReadStates()=0;
	virtual bool FilterIntensity(int I, int F)=0;

	bool degeneracy_filter(int gammaI,int gammaF);
	const char* branch(int jF,int jI);
	int GetNumberStates(){return eigenvalues.size();};
	int GetRecord(int state,int deg){return eigenvalues.at(state).irec.at(deg);};
	int GetRecordLength(int state){return eigenvalues.at(state).rec_len;};
	int GetRoot(int state,int deg){return eigenvalues.at(state).iroot.at(deg);};
	double GetEnergy(int state){return eigenvalues.at(state).energy;};
	int GetJ(int state){return eigenvalues.at(state).jval;};
	int GetJIndex(int state){return eigenvalues.at(state).jind;};
	int GetGamma(int state){return eigenvalues.at(state).igamma;};
	int GetLevel(int state){return eigenvalues.at(state).ilevel;};
	int GetNdeg(int state){return eigenvalues.at(state).ndeg;};
	void GetTransitionDetails(int & num_initial, int & num_trans,int & max_trans);


};
