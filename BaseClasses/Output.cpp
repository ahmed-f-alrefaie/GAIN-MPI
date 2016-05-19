#include "Output.h"
#include <cmath>
Output::Output(States* pstates,std::vector<double> p_gns,double ptemperature,double pZPE,double ppartition,double thresh,int pmaxDeg,bool red,bool do_file,const char* pfilename,bool full_line,bool compute_intens) : BaseProcess(){
	Log("Initializing output class\n");
	
	states = pstates;
	temperature = ptemperature;
	Q = ppartition;
	maxDeg=pmaxDeg;
	gns = p_gns;
	threshold=thresh;
	to_file = do_file;
	ZPE = pZPE;
	reduced = red;
	output_linestrengths = full_line;
	compute_intensities = compute_intens;
	
	if(pfilename != NULL && do_file){
		sprintf(filename,"%s__%d__.out",pfilename,m_process_id);
		to_file = true;
		
	}else{
		//printf("No filename!\n");
		to_file = false;
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

	char buffer[1024];
	char lsbuf[300];
	lsbuf[0]='\0';

	
	double energyI = states->GetEnergy(iLevelI);
	double nu_if = states->GetEnergy(iLevelF) - states->GetEnergy(iLevelI);
	energyI-=ZPE;
	if(nu_if < 1e-40)
		return;
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
	ls =0.0;
	for(int idegF=0; idegF < ndegF; idegF++){
		for(int idegI=0; idegI < ndegI; idegI++){
			if(!states->DegeneracyFilter(gammaI,gammaF,idegI,idegF))
				continue;
			linestr = linestrength[idegI + idegF*maxDeg];
			ls +=(linestr*linestr);
			if(output_linestrengths)
				sprintf(lsbuf + strlen(lsbuf)," %16.9E",linestr);		
			
		}
	}
	ls/=double(ndegI);
	if (reduced && ndegF!=1 && ndegI !=1 ) 
		ls*=double(ndegI);

	A_einst = ACOEF*double(2*jI+1)*ls*fabs(nu_if)*fabs(nu_if)*fabs(nu_if);

	if(A_einst < threshold)
		return;


	sprintf(buffer,"%12.6f %8d %4d %4d <- %8d %4d %4d %16.9E ",nu_if,indexF+1,jF,gammaF+1,indexI+1,jI,gammaI+1,A_einst);
	if(compute_intensities){
		double intens = CMCOEF*A_einst*gns[gammaF]*(2.0*double(jF) + 1.0)
							*std::exp(-beta*energyI)*(1.0-std::exp(-beta*nu_if))/(nu_if*nu_if*Q);
		sprintf(buffer + strlen(buffer), "%16.9E ",intens);
		
	}
	if(output_linestrengths){
		sprintf(buffer + strlen(buffer), "[ %s ] ",lsbuf);
	}


	fprintf(output,"%s || \n",buffer);
}

void Output::Close(){
	if(output != stdout && output != NULL){
		fclose(output);
	}
}
