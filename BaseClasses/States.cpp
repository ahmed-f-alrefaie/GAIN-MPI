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
	erange_upper =input.GetErangeUpper();
	ZPE = input.GetZPE();
	freq_window = input.GetFreqWindow();
	isym_pairs = input.GetISymPairs();
	igamma_pair = input.GetIGammaPair();
	reduced = input.IsReduced();
	nsize_max = 0;

}


bool States::filter_state(double energy,int igamma){

	return filter_lower(energy,igamma) || filter_upper(energy,igamma);

}
bool States::filter_lower(double energy,int igamma){
	return (energy - ZPE) <= erange_lower[1] && (energy - ZPE) >= erange_lower[0] && isym_do[igamma];
}
bool States::filter_upper(double energy,int igamma){
	return (energy - ZPE) <= erange_upper[1] && (energy - ZPE) >= erange_upper[0] && isym_do[igamma];


}
bool States::FilterIntensity(int I, int F){
	const EigenStates* eig_i = &eigenvalues.at(I);
	const EigenStates* eig_f = &eigenvalues.at(F);
	int gammaI = eig_i->igamma;
	int gammaF = eig_f->igamma;
	double nu_if = eig_f->energy - eig_i->energy;
//	Log("%12.6f %12.6f %d %d %12.6f <-  %12.6f %d %d %12.6f \n",nu_if,eig_f->energy,eig_f->jval,gammaF,gns[gammaF],eig_i->energy,eig_i->jval,gammaI,gns[gammaI]);
	int jI = eig_i->jval;
	int jF = eig_f->jval;
	return  gns[gammaI]> 0.0 &&
		gns[gammaF] > 0.0 &&
		nu_if >= freq_window[0] &&
		nu_if <= freq_window[1] &&
		filter_lower(eig_i->energy,gammaI) &&
		filter_upper(eig_f->energy,gammaF) &&
		isym_pairs[gammaI]==isym_pairs[gammaF] &&
		igamma_pair[gammaI]==gammaF &&
		abs(jF-jI)<=1		    &&
		jI+jF >=1;
		

}
bool States::FilterAnyTransitionsFromJ(int I, int J){

	int jI = eigenvalues.at(I).jval;
	if(abs(jI-J)>1 || jI+J==0)
		return false;
	for(int ilevelF=0; ilevelF < Neigenlevels; ilevelF++){
		if(eigenvalues.at(ilevelF).jval != J)
			continue;

		if(FilterIntensity(I, ilevelF))
			return true;

	}
	return false;

}


bool States::DegeneracyFilter(int gammaI,int gammaF,int idegI,int idegF){
	if(sym_degen[gammaI] == 1 || sym_degen[gammaF] == 1)
		return true;
	if(!reduced)
		return true;
	
	if(idegI==1 && idegF==0)
		return true;
	else
		return false;
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
	for(int ilevelI=0; ilevelI < Neigenlevels; ilevelI++){

		if(!FilterLowerState(ilevelI)) continue;
		per_initial_states = 0;
		for(int ilevelF=0; ilevelF < Neigenlevels; ilevelF++){

			if(!FilterIntensity(ilevelI,ilevelF)) continue;

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


}
