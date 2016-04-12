#include "Output.h"

Output::Output(States* pstates,double ptemperature,double ppartition,double thresh,int pmaxDeg,const char* pfilename) : BaseProcess(){
	Log("Initializing output class\n");
	
	states = pstates;
	temperature = ptemperature;
	Q = ppartition;
	maxDeg=pmaxDeg;
	threshold=thresh;
	to_file = false;

	if(pfilename != NULL){
		sprintf(filename,"%s__%d__.out",pfilename,m_process_id);
		to_file = true;
		
	}

	beta = PLANCK * VELLGT / (BOLTZ * temperature);

}
	
void Output::Initialize(){
	if(to_file){
		VerboseLog("Opeing output file: %s\n",filename);
		 output = fopen(filename,"wb");
		
	}else{
		output = stdout;
	}
	if(output == NULL){
		LogErrorAndAbort("Problem opening either stdout or file");

	}


	


}

void Output::OutputLinestrength(int iLevelI,int iLevelF,double* linestrength){
	double ls = 0.0;
	double linestr = 0.0;
	double A_einst = 0.0;

	double nu_if = states->GetEnergy(iLevelF) - states->GetEnergy(iLevelI);
	//Lower
	int jI= states->GetJ(iLevelI);
	int gammaI = states->GetGamma(iLevelI);
	int ndegI = states->GetNdeg(iLevelI); 
	int indexI = states->GetLevel(iLevelI);
	//Upper
	int jF= states->GetJ(iLevelF);
	int gammaF = states->GetGamma(iLevelF);
	int ndegF = states->GetNdeg(iLevelF); 
	int indexF = states->GetLevel(iLevelF);

	for(int idegF=0; idegF < ndegF; idegF++){
		for(int idegI=0; idegI < ndegI; idegI++){
			linestr = linestrength[idegI + idegF*maxDeg];
					
			ls +=(linestr*linestr);
		}
	}

			

	ls/=double(ndegI);

	A_einst = ACOEF*double((2*jI)+1)*ls*nu_if*nu_if*abs(nu_if);

	if(A_einst < threshold)
		return;

	fprintf(output,"%12.6f %8d %4d %4d <- %8d %4d %4d %16.8E \n",nu_if,indexF+1,jF,gammaF+1,indexI+1,jI,gammaI+1,A_einst);
}

void Output::Close(){
	if(output != stdout && output != NULL){
		fclose(output);
	}
}
