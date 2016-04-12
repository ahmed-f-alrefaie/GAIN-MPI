
#include "TroveInput.h"

void TroveInput::ReadInput(const char* filename){

	std::string line;
	stream.open(filename);
	if(!stream.good()){
		Log("Could not open file %s\n",filename);
		exit(0);
	}
	std::string clean_line;
	while(std::getline(stream,line)){
		//Output the line		
		clean_line = PrepareString(line);
		
		Log("%s\n",clean_line.c_str());
		//Split the line
		std::vector<std::string> split_line = split(clean_line);
				
		if(split_line.size()==0)
			continue;
		
		std::string input_flag = trim(split_line[0]);
		
		if(input_flag == "MEM"){
			memory = atoi(split_line[1].c_str());
			if(split_line.size()==3){
				if(trim(split_line[2])=="GB")
					memory*=GIGABYTE_TO_BYTE;
				else if(trim(split_line[2])=="MB")
					memory*=MEGABYTE_TO_BYTE;
			}
			
			BaseManager::InitializeGlobalMemory(memory);
			
		}else if (input_flag == "SYMGROUP"){
			symmetry_group = trim(split_line[1]);
		}else if(input_flag == "INTENSITY"){
			ReadIntensity();
		}
		
		
		
	}
	maxJ=0;
	HandleSymmetery();
	ComputeIGammaPair();
	erange.push_back(std::min(erange_lower[0],erange_upper[0]));
	erange.push_back(std::max(erange_lower[1],erange_upper[1]));

	for(int i = 0; i < jVals.size(); i++)
		maxJ = std::max(maxJ,jVals[i]);

	Log("Energy Range %12.6f - %12.6f\n  Jmax: %d\n",erange[0],erange[1],maxJ);
	Log("Lower Range %12.6f - %12.6f\n",erange_lower[0],erange_lower[1]);
	Log("Upper Range %12.6f - %12.6f\n",erange_upper[0],erange_upper[1]);
	Log("total Memory %12.6f GB",double(memory)/double(GIGABYTE_TO_BYTE));

}

std::string TroveInput::PrepareString(std::string line){
	char linestr[1024];
		//Output the line		
		//Copy
	strcpy (linestr, line.c_str());
	for(int i = 0; i < line.length(); i++)
		linestr[i]=toupper(linestr[i]);
		//Remove all of the brackets
	RemoveBrackets(linestr);
	return std::string(linestr);
}

void TroveInput::HandleSymmetery(){
	sym_maxdegen = 1;
	Log("SYMMETRY-STUFF");
	std::cout<<symmetry_group<<std::endl;
	if(symmetry_group=="C2V"){
		
		sym_degen.push_back(1);

		sym_degen.push_back(1);

		sym_degen.push_back(1);

		sym_degen.push_back(1);

	}else if(symmetry_group=="C3V"){
		sym_degen.push_back(1);

		sym_degen.push_back(1);
		
		sym_degen.push_back(2);
	}else if(symmetry_group=="D2H"){
		sym_degen.push_back(1);


		sym_degen.push_back(1);

		
		sym_degen.push_back(1);


		sym_degen.push_back(1);


		sym_degen.push_back(1);

		
		sym_degen.push_back(1);


		sym_degen.push_back(1);

		
		sym_degen.push_back(1);

	}
	

	sym_nrepres = sym_degen.size();
	for(int i = 0; i < sym_nrepres; i++)
		sym_maxdegen = std::max(sym_maxdegen,sym_degen[i]);
	Log("Max degeneracy = %d\n",sym_maxdegen);
}

void TroveInput::RemoveBrackets(char* string){
	int len =strlen(string);
	for(int i = 0; i < len; i++){
		if(string[i]==',')
			string[i]=' ';
		if(string[i]=='('){
			for(int j = i; j < len; j++){
				if(string[j]==')'){
					string[j]=' ';
					break;
				}
				string[j]=' ';
				
			}
		}
			
	}

}

void TroveInput::ComputeIGammaPair(){

	int ngamma=0;
	//igamma_pair.resize(
	Log("Compute Gamma pairs\n");
	for(int igammaI = 0; igammaI < sym_nrepres; igammaI++)
	{
		ngamma=0;
		igamma_pair.push_back(igammaI);
		Log("igammaI = %i sym_nrepres = %i\n",igammaI,sym_nrepres);
		for(int igammaF = 0; igammaF < sym_nrepres; igammaF++)
		{
			if(igammaI!=igammaF && isym_pairs[igammaI]==isym_pairs[igammaF])
			{
				igamma_pair[igammaI] = igammaF;
				ngamma++;
				if(ngamma>1){
					LogErrorAndAbort("find_igamma_pair: Assumption that selection rules come in pairs is not fulfilled!\n");
					
					
				}
			}
				
		}
		//
		if(gns[igammaI] != gns[igamma_pair[igammaI]]){
			LogErrorAndAbort("find_igamma_pair: selection rules do not agree with Gns\n");
					
		}
		
	}



}


void TroveInput::ReadIntensity(){

	std::string line;
	while(std::getline(stream,line)){
		
		std::string clean_line = PrepareString(line);
		Log("%s\n",clean_line.c_str());
		std::vector<std::string> split_line = split(clean_line);
		
		std::string input_flag = trim(split_line[0]);
		if(input_flag == "END")
			return;

		if(input_flag == "TEMPERATURE")
			temperature = atof(split_line[1].c_str());
		else if(input_flag == "QSTAT")
			q_stat = atof(split_line[1].c_str());
		else if(input_flag == "ZPE")
			ZPE = atof(split_line[1].c_str());
		else if(input_flag == "J"){
			int start_J = atoi(split_line[1].c_str());
			int end_J = atoi(split_line[2].c_str());
			nJ = abs(end_J-start_J + 1);
			if(nJ < 2 || start_J > end_J){
				Log("cannot do intensities with configuration J %i,%i\n",start_J,end_J);
				exit(0);
			}
			for(int i = 0; i < nJ; i++){
				jVals.push_back(start_J+i);
			}
			
		}else if(input_flag == "FREQ-WINDOW"){
			freq_window.push_back(atof(split_line[1].c_str()));
			freq_window.push_back(atof(split_line[2].c_str()));
		}else if(input_flag == "ENERGY"){
			//Log("Reading ENERGY");
			for(int i = 1; i < split_line.size(); i++){
				if(trim(split_line[i])=="LOW" || trim(split_line[i])=="LOWER"){
					erange_lower.push_back(atof(split_line[++i].c_str()));
					erange_lower.push_back(atof(split_line[++i].c_str()));
				}
				if(trim(split_line[i])=="UPPER"){
					erange_upper.push_back(atof(split_line[++i].c_str()));
					erange_upper.push_back(atof(split_line[++i].c_str()));
				}



			}		

		
		}else if(input_flag == "GNS"){
			for(int i = 1; i < split_line.size(); i++){
				double t_gns = atof(split_line[i].c_str());
				gns.push_back(atof(split_line[i].c_str()));
				if(t_gns< 1e-15){
					isym_do.push_back(false);
				}else{
					isym_do.push_back(true);
				}
		
			}
			
		}else if(input_flag == "SELECTION"){
			for(int i = 1; i < split_line.size(); i++){
				isym_pairs.push_back(atoi(split_line[i].c_str()));
			}			
		}else if(input_flag == "THRESH_LINE"){
			thresh_linestrength = atof(split_line[1].c_str()); 
		}else if ( input_flag == "SYMMETRY"){
			if(trim(split_line[1])=="REDUCED"){
				reduced=true;
			}
		}




	}


}

