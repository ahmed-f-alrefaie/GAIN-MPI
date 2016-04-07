#include "States.h"





States::States(Input & input): BaseProcess(){
	
	Neigenlevels=0;
	Neigenroots=0;

	sym_maxdegen = input.GetSymMaxDegen();
	minJ = 100000;
	maxJ=0;

	jVals = input.GetJvals();
	//sym_gamma = input.GetSymGamma();
	for(int i = 0; i < jVals.size(); i++){
		minJ = std::min(jVals[i],minJ);
		maxJ = std::max(jVals[i],maxJ);
	}

	gns = input.GetGNS();
	isym_do=input.GetISymDo();
	sym_degen=input.GetSymmetryDegen();
	sym_nrepres = input.GetNSym();
	erange=input.GetErange();
	erange_lower =input.GetErangeLower();
	erange_upper =input.GetErangeLower();
	ZPE = input.GetZPE();
	freq_window = input.GetFreqWindow();
	isym_pairs = input.GetISymPairs();
	igamma_pair = input.GetIGammaPair();

	nsize_max = 0;

}


bool States::filter_state(double energy,int igamma){

	return filter_lower(energy,igamma) || filter_upper(energy,igamma);

}
bool States::filter_lower(double energy,int igamma){
	return (energy - ZPE) <= erange_lower[1] && (energy - ZPE) >= erange_lower[0] && isym_do[igamma] > 0.0;
}
bool States::filter_upper(double energy,int igamma){
	return (energy - ZPE) <= erange_upper[1] && (energy - ZPE) >= erange_upper[0] && isym_do[igamma] > 0.0;


}



/*bool States::degeneracy_filter(int gammaI,int gammaF){
	if(sym_degen[gammaI]==1 || sym_degen[gammaF]==1) return true;
 

}
*/

void States::GetTransitionDetails(int & num_initial, int & num_trans,int & max_trans){
	int no_initial_states=0;
	int transitions = 0;
	int max_trans_per_initial = 0;
	
	int per_initial_states = 0;
	Log("\n\nPredicting how many transitions to compute\n");
/*	for(int ilevelI=0; ilevelI < Neigenlevels; ilevelI++){

		if(!energy_filter_lower(ilevelI)) continue;
		per_initial_states = 0;
		for(int ilevelF=0; ilevelF < Neigenlevels; ilevelF++){

			if(!energy_filter_upper(ilevelF)) continue;

			if(!intensity_filter(ilevelI,ilevelF)) continue;

			transitions++;
			per_initial_states++;

		}
		max_trans_per_initial = std::max(max_trans_per_initial,per_initial_states);
		no_initial_states++;


	}
	Log("\n.....We have counted %i transitions with %i initial states and largest no transitions is %i\n\n",int(transitions),no_initial_states,max_trans_per_initial);
	
	num_initial = no_initial_states;
	num_trans = transitions;
	max_trans = max_trans_per_initial;
*/

}
