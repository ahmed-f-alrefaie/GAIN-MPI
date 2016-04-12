#include "TroveBasisSet.h"



TroveBasisSet* TroveBasisSet::j0BasisSet = NULL;






TroveBasisSet::TroveBasisSet(int J,int sym_n, std::vector<int> psym_degen) : BasisSet(J){ 
	sym_nrepres = sym_n;
	sym_degen=psym_degen;

}


void TroveBasisSet::Initialize(){
	ReadBasisSet();
	Correlate();
	ComputeijTerms();
}

void TroveBasisSet::ReadBasisSet(){
	//int* Ntotal;
	char filename[1024];
	std::string line;
	char* line_ptr;
	int ncases,nlambdas,ncontr,nclasses,icontr,iroot;

	//Null pointer to be safe :3c
	icontr2icase=NULL;
    	icase2icontr=NULL;
   	ktau=NULL;
    	k=NULL;
    	contractive_space=NULL;
    	index_deg=NULL;
    	rot_index=NULL;
    	icontr_correlat_j0=NULL;
    	iroot_correlat_j0=NULL;
    	nsize=NULL;
    	irr=NULL;
    	Ntotal=NULL;
	
	
	
	//Get the filename
	sprintf(filename,"j0eigen_quanta%d.chk",jval);
	
	Log("Reading %s\n",filename);

	//Open file
	std::ifstream eig_qu(filename);
	
	//Begin reading
	std::getline(eig_qu,line);
	
	if(trim(line).compare("Start Primitive basis set")!=0)
	{
		LogErrorAndAbort("[bset_contr_factory]: bad header %s",filename);
		//fLog(stderr,"[bset_contr_factory]: bad header");
	}	

	std::getline(eig_qu,line);
	/*
	       read(iounit, '(4i8)') ncases, nlambdas, ncontr, nclasses

       bset_contr(jind)%Maxsymcoeffs = ncases
       bset_contr(jind)%max_deg_size = nlambdas
       bset_contr(jind)%Maxcontracts = ncontr
       bset_contr(jind)%nclasses     = nclasses
	*/
	//jval = J;
	Maxsymcoeffs  = ncases = strtol(line.c_str(),&line_ptr,0);
	max_deg_size = nlambdas = strtol(line_ptr,&line_ptr,0);
	Maxcontracts = ncontr = strtol(line_ptr,&line_ptr,0);
	Nclasses = nclasses = strtol(line_ptr,&line_ptr,0);

	Maxcontracts = ncontr;
	//Log("MAx_deg_size = %d\n",max_deg_size);
	
	dimenMax = std::max(dimenMax,Maxcontracts);
	//Allocate memory
	/*
	       allocate(bset_contr(jind)%index_deg(ncases),bset_contr(jind)%contractive_space(0:nclasses, ncases),stat = info)
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%contractive_space),kind(bset_contr(jind)%contractive_space))
       !
       allocate(bset_contr(jind)%nsize(sym%Nrepresen),stat = info)
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%nsize),kind(bset_contr(jind)%nsize))
       !
       allocate(bset_contr(jind)%icontr2icase(ncontr, 2),bset_contr(jind)%icase2icontr(ncases, nlambdas),stat = info)
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%icontr2icase),kind(bset_contr(jind)%icontr2icase))
       call ArrayStart('bset_contr',info,size(bset_contr(jind)%icase2icontr),kind(bset_contr(jind)%icase2icontr))
       */
       index_deg = new TO_PTintcoeffsT[ncases];
       contractive_space = new int[(nclasses+1)*ncases];
       nsize = new int[sym_nrepres];
       icontr2icase = new int[ncontr*2];
       icase2icontr = new int[ncases*nlambdas];
	//cout<<nsize[0]<<endl;
       	#ifdef DEBUG
		Log("%i %i %i %i\n",Maxsymcoeffs,max_deg_size,Maxcontracts,Nclasses);
	#endif
       icontr = 0;
       iroot = 0;
       
       for(int icase = 0; icase < ncases; icase++)
       {
       		getline(eig_qu,line);
       		//read(iounit,*) nlambdas, bset_contr(jind)%contractive_space(0:nclasses, icase)
        	//bset_contr(jind)%index_deg(icase)%size1 = nlambdas
        	nlambdas = strtol(line.c_str(),&line_ptr,0);
        	for(int i = 0; i <=nclasses; i++)
        		 contractive_space[i + icase*(nclasses+1)] = strtol(line_ptr,&line_ptr,0)-1;
        	
        	index_deg[icase].size1 = nlambdas;
		#ifdef DEBUG
			Log("nlambdas = %i contr_space_1 = %i contr_space_2 = %i \n", nlambdas, contractive_space[0 + icase*(nclasses+1)], contractive_space[1 + icase*(nclasses+1)]);
		#endif
        	
        	//    allocate(bset_contr(jind)%index_deg(icase)%icoeffs(0:nclasses, nlambdas))
        	index_deg[icase].icoeffs = new int[(nclasses+1)*nlambdas]; 
        	
        	for(int ilambda = 0; ilambda < nlambdas; ilambda++)
        	{
        		std::getline(eig_qu,line);
        		//read(iounit, '(i8, i6, <nclasses>i6)') iroot, bset_contr(jind)%index_deg(icase)%icoeffs(0:nclasses, ilambda)
        		iroot = strtol(line.c_str(),&line_ptr,0)-1;
        		for(int i = 0; i <=nclasses; i++)
        			index_deg[icase].icoeffs[i+ ilambda*(nclasses+1)] = strtol(line_ptr,&line_ptr,0)-1;
        	
        		#ifdef DEBUG
				Log("iroot = %i icoeffs_1 = %i icoeffs_2 = %i \n", iroot, index_deg[icase].icoeffs[0+ ilambda*(nclasses+1)], index_deg[icase].icoeffs[1+ ilambda*(nclasses+1)]);
			#endif	
        		

             		icontr2icase[icontr + 0*ncontr]      = icase;
             		icontr2icase[icontr + 1*ncontr]      = ilambda;
        		#ifdef DEBUG
				Log("icontr = %i icase = %i ilambda = %i \n", icontr,icase,ilambda);
			#endif	
             		icase2icontr[icase + ilambda*ncases] = icontr;
        	
        		
        		
        		if(iroot != icontr)
        		{
        			LogErrorAndAbort("[bset_contr_factory] wrong indexing icontr = %i iroot = %i\n",icontr,iroot);
        			
        		}
        		icontr = icontr + 1;
        	}
        	  
       }
	//if(jval==1){Log("icontr2icase[:,2]\n");     
	//for(int i = 0; i < ncontr; i++)
	//	Log("%i\n",icontr2icase[i + 1*ncontr]);

	//exit(0);
		
	//}

	if (icontr != ncontr)
	{
		LogErrorAndAbort("[bset_contr_factory] wrong indexing\n");
		
	}
	
	//read(iounit, '(2i8)') ncases, nlambdas
	std::getline(eig_qu,line);
	ncases = ncases = strtol(line.c_str(),&line_ptr,0);
	nlambdas =nlambdas = strtol(line_ptr,&line_ptr,0);
	
	//allocate(bset_contr(jind)%rot_index(ncases, nlambdas), stat = info) 
	rot_index=new TO_PTrotquantaT[ncases*nlambdas];
	Log("j= %i rot_index size =%i,ncases=%i, nlambdas = %i\n",jval,ncases*nlambdas,ncases,nlambdas); 
	int icase,ilambda;
	while(std::getline(eig_qu,line))
	{
		if( trim(line).compare("End Primitive basis set")==0)
			break;
		/*read(buf, *) icase, ilambda, bset_contr(jind)%rot_index(icase, ilambda)%j,                                       &
                                       bset_contr(jind)%rot_index(icase, ilambda)%k,                                       &
                                       bset_contr(jind)%rot_index(icase, ilambda)%tau
		
		*/
		icase = strtol(line.c_str(),&line_ptr,0)-1;
		ilambda = strtol(line_ptr,&line_ptr,0)-1;
		rot_index[icase + ilambda*ncases].j = strtol(line_ptr,&line_ptr,0);
		rot_index[icase + ilambda*ncases].k = strtol(line_ptr,&line_ptr,0);
		rot_index[icase + ilambda*ncases].tau = strtol(line_ptr,&line_ptr,0);
			
		#ifdef DEBUG
			Log("icase = %i ilambda = %i j = %i k = %i tau = %i\n",icase,ilambda,rot_index[icase + ilambda*ncases].j,rot_index[icase + ilambda*ncases].k,rot_index[icase + ilambda*ncases].tau);
		#endif
	
	}
	
	
		//Begin reading
	std::getline(eig_qu,line);
	//cout<<line<<endl;
	if(trim(line).compare("Start irreducible transformation")!=0)
	{
		LogErrorAndAbort("[bset_contr_factory]: wrong sym-footer");
		
	}
	
	irr = new TO_PTrepresT[sym_nrepres];
	Ntotal = new int[sym_nrepres];
	getline(eig_qu,line);

	mat_size = strtol(line.c_str(),&line_ptr,0);

	Log("matsize = %i\n",mat_size);

	//cout<<line<<endl;
	//mat_size = mat_size;

	//read(iounit,*) Ntotal(1:sym%Nrepresen)
	std::getline(eig_qu,line);

	Ntotal[0] = strtol(line.c_str(),&line_ptr,0); 

	//Log("Ntotal[0]=%i\n",Ntotal[0]);

	//cout<<line<<endl;
	for(int i = 1; i < sym_nrepres; i++){

		Ntotal[i] = strtol(line_ptr,&line_ptr,0);
		//Log("Ntotal[%i]=%i\n",i,Ntotal[i]);
	}
	//(bset_contr(jind)%irr(igamma)%N(bset_contr(jind)%Maxsymcoeffs),bset_contr(jind)%irr(igamma)%repres(Ntotal(igamma),sym%degen(igamma),mat_size)
	//Log("nrepres = %i",sym_nrepres);
	for(int igamma = 0; igamma < sym_nrepres; igamma++)
	{
		irr[igamma].N = new int[Maxsymcoeffs];
		irr[igamma].repres = new double[Ntotal[igamma]*sym_degen[igamma]*mat_size];
		//Log("J=%d gamma=%d Ntotal=%d Sym_degen=%d mat_size=%d repres_size=%d\n",
		 //   jval,   igamma,    Ntotal[igamma],    sym_degen[igamma],mat_size,   (Ntotal[igamma]*sym_degen[igamma]*mat_size));
		
		//           do icoeff = 1,bset_contr(jind)%Maxsymcoeffs
                //		!
                // 		read(iounit,*) bset_contr(jind)%irr(igamma)%N(icoeff)
                //		!
                //	     enddo
                
                for(int icoeff = 0; icoeff <  Maxsymcoeffs; icoeff++)
                {
                		std::getline(eig_qu,line);
				//std::cout<<line<<std::endl;
				irr[igamma].N[icoeff] = strtol(line.c_str(),&line_ptr,0);
		}
                

              for(int icoeff = 0; icoeff <  Ntotal[igamma]; icoeff++)
              {
              		for(int ideg =0; ideg <  sym_degen[igamma]; ideg++)
              		{
              			std::getline(eig_qu,line);	
				//std::cout<<line<<std::endl;
              			irr[igamma].repres[icoeff + ideg*Ntotal[igamma]] = strtod(line.c_str(),&line_ptr);
				#ifdef DEBUG
					Log("irr[%i].repres[%i,%i,%i]=%11.4e\n",igamma,icoeff,ideg,0,irr[igamma].repres[icoeff + ideg*Ntotal[igamma] + 0*Ntotal[igamma]*sym_degen[igamma]]);
				#endif
              			for(int i = 1; i < mat_size; i++){
					irr[igamma].repres[icoeff + ideg*Ntotal[igamma] + i*Ntotal[igamma]*sym_degen[igamma]] = strtod(line_ptr,&line_ptr);
					#ifdef DEBUG
						Log("irr[%i].repres[%i,%i,%i]=%11.4e\n",igamma,icoeff,ideg,i,irr[igamma].repres[icoeff + ideg*Ntotal[igamma] + i*Ntotal[igamma]*sym_degen[igamma]]);
					#endif
				}
			}
	      }

		
			
	}
	std::getline(eig_qu,line);
	//std::cout<<line<<std::endl;
	if(trim(line).compare("End irreducible transformation")!=0)
	{
		LogErrorAndAbort("[bset_contr_factory]: wrong irrep-footer");
		
	}
	
	Log("Done!\n");

	
}

void TroveBasisSet::Correlate(){

	if(j0BasisSet == NULL){
		j0BasisSet = new TroveBasisSet(0,sym_nrepres,sym_degen);
		j0BasisSet->Initialize();

	}


	 Log("Establish the correlation between the indexes of J=0 and J=%i contr. basis funct.\n",jval);	
	
	if(j0BasisSet->jval != 0)
	{
		 LogErrorAndAbort("[correlate_index] index_correlation: bset_contrj0 is not for J=0\n");
		
	}
	
	int nclasses = Nclasses;
	
	if(j0BasisSet->Nclasses != nclasses)
	{
		//fprintf(stderr,"[index_correlation]: Nclasses are different for diff. J");
		//exit(0);
	}
	
	int icase,jcase,ilambda,jlambda,icontr,jcontr,info, iroot,jroot;
	int ilevel,ideg,t_k,t_tau,dimen,irow,icol;

	int*	cnu_i, *cnu_j;

	bool found;
	
	//allocate(cnu_i(1:nclasses),cnu_j(1:nclasses),stat = info)
	cnu_i = new int[nclasses];
	cnu_j = new int[nclasses];
	 Log("Maxsymj0 = %i %i %i \n",j0BasisSet->Maxcontracts,j0BasisSet->Maxsymcoeffs,	j0BasisSet->max_deg_size);
	 Log("Maxsym = %i %i %i \n",Maxcontracts,Maxsymcoeffs,max_deg_size);
	 Log("J0 from this J: = %d\n",Maxcontracts/(2*jval + 1));
	//allocate(bset_contr(jind)%icontr_correlat_j0(bset_contr(jind)%Maxsymcoeffs,bset_contr(jind)%max_deg_size), stat = info)
	icontr_correlat_j0 = new int[Maxsymcoeffs*max_deg_size];
	// Log("Do I root correlate\n");
	// Log("Start_correlation\n");
	//bool contract_space;
	//bool index_deg;
	//do icase = 1, bset_contr(jind)%Maxsymcoeffs
	for(icase = 0; icase < Maxsymcoeffs; icase++)
	{	
 		// Log("icase=%i\n",icase);
		//cnu_i(1:nclasses) = bset_contr(jind)%contractive_space(1:nclasses, icase)
		for(int c_i = 0; c_i < nclasses; c_i++)
			cnu_i[c_i] = contractive_space[(c_i+1) + icase*(nclasses+1)];
		//memcpy(cnu_i,contractive_space + icase*(nclasses+1) + 1,sizeof(int)*nclasses);
		//do ilambda = 1, bset_contr(jind)%index_deg(icase)%size1
		for(ilambda = 0; ilambda < index_deg[icase].size1; ilambda++)
		{
			found = false;

			//do jcase = 1, bset_contr(1)%Maxsymcoeffs
			for(jcase = 0; jcase < j0BasisSet->Maxsymcoeffs; jcase++)
			{
				//cnu_j(1:nclasses) = bset_contr(1)%contractive_space(1:nclasses, jcase)
				for(int c_j = 0; c_j < nclasses; c_j++)
					cnu_j[c_j] = j0BasisSet->contractive_space[(c_j+1) + jcase*(nclasses+1)];
				//memcpy(cnu_j,j0BasisSet->contractive_space + (1 + jcase*(nclasses+1)),sizeof(int)*nclasses);
				//do jlambda = 1, bset_contr(1)%index_deg(jcase)%size1
				for(jlambda = 0; jlambda < j0BasisSet->index_deg[jcase].size1; jlambda++)
				{
					/*			                 if (all(cnu_i(:) == cnu_j(:))  .and.   &
                   				  all(bset_contr(   1)%index_deg(jcase)%icoeffs(1:nclasses,jlambda) == &
                  				       bset_contr(jind)%index_deg(icase)%icoeffs(1:nclasses,ilambda))) then 
                  				       !
                  				       found = .true.
                 			        exit l_jcase
                			         !
                 					endif
               					  */  
					//for(int i = 0; i < nclasses; i++)
					//	 Log(" %i ",cnu_i[i]);
					// Log("\t|\t");
					//for(int i = 0; i < nclasses; i++)
					//	 Log(" %i ",cnu_j[i]);
					// Log("\n");					
					if(memcmp(cnu_i,cnu_j,sizeof(int)*nclasses)==0 && memcmp(j0BasisSet->index_deg[jcase].icoeffs + jlambda*(nclasses+1) + 1,
												 index_deg[icase].icoeffs + ilambda*(nclasses+1) + 1,
												 sizeof(int)*nclasses)==0)
					{
						found = true;
						break;
					}
				}//jlambda
				if(found) break;
			}//jcase
			if(!found)
			{
				 Log("icase =%i ilambda = %i\n",icase,ilambda);
				 LogErrorAndAbort("[index_correlation] not found for J = %i -> problems with checkpoints?\n", jval);
				//f Log(stderr,"[index_correlation] No correlation for J = %i\n", jval);
				
			}
			// Log("Done!");
			//jcontr = bset_contr(1)%icase2icontr(jcase,jlambda)
			jcontr = j0BasisSet->icase2icontr[jcase + jlambda*j0BasisSet->Maxsymcoeffs];
			// Log("jcontr = %i\n",jcontr);
			icontr_correlat_j0[icase + ilambda*Maxsymcoeffs] = jcontr;
			//             jcontr = bset_contr(1)%icase2icontr(jcase,jlambda)
            // !
            // bset_contr(jind)%icontr_correlat_j0(icase, ilambda) = jcontr
		}//ilambda
	}//icase
	// Log("Finished correlation\n");
//	       allocate(bset_contr(jind)%iroot_correlat_j0(bset_contr(jind)%Maxcontracts), stat = info)
//       call ArrayStart('bset_contr',info,size(bset_contr(jind)%iroot_correlat_j0),kind(bset_contr(jind)%iroot_correlat_j0))
//       allocate(bset_contr(jind)%ktau(bset_contr(jind)%Maxcontracts), stat = info)
//       call ArrayStart('bset_contr',info,size(bset_contr(jind)%ktau),kind(bset_contr(jind)%ktau))
//       allocate(bset_contr(jind)%k(bset_contr(jind)%Maxcontracts), stat = info)
//       call ArrayStart('bset_contr',info,size(bset_contr(jind)%k),kind(bset_contr(jind)%k))
	// Log("Allocate iroot j0\n");
	iroot_correlat_j0 = new int[Maxcontracts];
	// Log("Allocate ktau\n");	
	ktau = new int[Maxcontracts];
	//tau = new int[Maxcontracts];
	// Log("Allocate k\n");
	k = new int[Maxcontracts];
  	k_block_size= new int[jval+1];
   	kstart= new int[jval+1];
	kstart[0]=0;
	
/*	       do iroot = 1, bset_contr(jind)%Maxcontracts
          !
          icase   = bset_contr(jind)%icontr2icase(iroot, 1)
          ilambda = bset_contr(jind)%icontr2icase(iroot, 2)
          !
          jcontr = bset_contr(jind)%icontr_correlat_j0(icase, ilambda)
          !
          bset_contr(jind)%iroot_correlat_j0(iroot) = jcontr
          !
          ilevel  = bset_contr(jind)%contractive_space(0,icase)
          ideg    = bset_contr(jind)%index_deg(icase)%icoeffs(0,ilambda)
          !
          k      = bset_contr(jind)%rot_index(ilevel,ideg)%k
          tau    = bset_contr(jind)%rot_index(ilevel,ideg)%tau
          !
          bset_contr(jind)%ktau(iroot) = 2*k+tau
          bset_contr(jind)%k(iroot)    = k
*/
	int ncontr = Maxcontracts;
	int last_k_b = 0;
	int ksize=0;
	// Log("Doing k j zero stuf ncases = %i nlambda=%i\n",ncases,max_deg_size);
	for(int iroot = 0; iroot < Maxcontracts; iroot++)
	{
		
		//if(jval==23&&iroot > 14346) Log("Reached beginning\n");
          	icase   = icontr2icase[iroot];
		//if(jval==23&&iroot > 14346) Log("icase =%i\n",icase);
         	ilambda = icontr2icase[iroot + ncontr];
		//if(jval==23&&iroot > 14346) Log("ilambda =%i\n",ilambda);			
		jcontr = icontr_correlat_j0[icase + ilambda*Maxsymcoeffs];
		//if(jval==23&&iroot > 14346) Log("jcontr =%i\n",jcontr);
		iroot_correlat_j0[iroot] = jcontr;
		//if(jval==23&&iroot > 14346) Log("done\n");

		ilevel  = contractive_space[0 + icase*(nclasses+1)];
		//if(jval==23&&iroot > 14346) Log("ilevel = %i\n",ilevel);
		ideg    = index_deg[icase].icoeffs[0 + (nclasses+1)*ilambda];
		// Log("ideg = %i\n",ideg);
		//if(jval==23&&iroot > 14346) Log("rot parts\n");
		t_k      = rot_index[ilevel+ideg*ncases].k;
		//if(jval==23&&iroot > 14346) Log("k = %i\n",k);
          	t_tau    = rot_index[ilevel+ideg*ncases].tau;
		//if(jval==23&&iroot > 14346) Log("tau=%i\n",tau);
		ktau[iroot] = 2*t_k+t_tau & 1;
		//tau[iroot] = ktau[iroot] & 1;
          	k[iroot]    = t_k;
		if(last_k_b != t_k){
			k_block_size[last_k_b]=ksize;
			kstart[t_k]=kstart[last_k_b]+ksize;
			ksize=0;
			last_k_b=t_k;
		}

		ksize++;
		if(iroot_correlat_j0[iroot] >= j0BasisSet->Maxcontracts){
			 LogErrorAndAbort("Error in root correlats %i %i \n",iroot_correlat_j0[iroot],j0BasisSet->Maxcontracts);
			
		}
          	//#ifdef DEBUG
          	//Log("iroot = %i ilevel = %i icase = %i ideg=%i ilambda = %i jcontr = %i k = %i, tau = %i\n",iroot,ilevel,icase,ideg,ilambda,jcontr,t_k,t_tau); fflush(0);
          	//#endif 
		// Log("iroot = %i ilevel = %i icase = %i ideg=%i ilambda = %i jcontr = %i k = %i, tau = %i\n",iroot,ilevel,icase,ideg,ilambda,jcontr,k,tau);
		// Log("iroot = %i j0root = %i\n",iroot,jcontr);
          	
         }
	k_block_size[last_k_b] = ksize;
	Log("---------K block analysis for J=%d------------\n",jval);
	for(int i =0; i < jval+1; i++){
		Log("K block %d start: %d end %d size: %d \n",i,kstart[i],kstart[i]+k_block_size[i]-1,k_block_size[i]);
	}
	

	// Log("Cleaning up\n");
	delete[] cnu_j;
	delete[] cnu_i;
	 Log("done\n");	


}
void TroveBasisSet::ComputeijTerms(){

    Log("Computing ijterms for J= %i\n",jval);
    ijterms = new int[Maxsymcoeffs*sym_nrepres];
    
    for(int igammaI= 0; igammaI < sym_nrepres; igammaI++)
    {
    	int Nterms = 0;
    	for(int iterm = 0; iterm < Maxsymcoeffs; iterm++)
    	{
    		#ifdef DEBUG
			printf("iterm,igamma = %i,%i Nterms = %i!\n",iterm,igammaI,Nterms);
    		#endif 
    		ijterms[iterm + igammaI*Maxsymcoeffs] = Nterms;
    		Nterms += irr[igammaI].N[iterm];
    	}
    }
    
    Log(".....Done!\n");


}

