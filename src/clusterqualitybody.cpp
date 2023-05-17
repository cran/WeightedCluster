#ifdef CLUSTERQUALITY_INCLUDED

#if !defined(DISTOBJECT_VERSION) && !defined(DISTMATRIX_VERSION)
	#error Clustequality version not defined
#endif

void CLUSTERQUALITY_FUNCNAME(double * distmatrix, int * clusterid, double *weights, int nelements, double* stats, int nclusters, double * errors2, KendallTree &kendall){
	TMRLOG(2,"Computing statitstics\n");
	double totweights=0, wxy=0, wx=0, wy=0, wx2=0, ww, xx, covxy, covx, covy, pearson, xb, yb, xw, xxw;
	//double totweights=0, wxy=0,wxy2=0, wx=0, wy=0, wx2=0, ww, xx, covxy, covx, covy, pearson, xb, yb, xw, xxw;
	int i, j, ij=0, iclustIndex;
	double *errors = (double*) R_alloc(nclusters,sizeof(double));
	int *tablesizes = (int*) R_alloc(nclusters,sizeof(int));
	double *sizes = (double*) R_alloc(nclusters,sizeof(double));
	for(i=0;i <nclusters;i++){
		tablesizes[i]=-1;
		errors2[i]=0;
		errors[i]=0;
		sizes[i]=0;
	}
	
	//Maybe use  R_qsort_I to sort distances between individuals first
	CmpCluster * cmpclust=NULL;
	CmpCluster * ZeroDist;
	KendallTreeIterator it;
	it=kendall.find(0.0);
	if (it!=kendall.end()) {
		ZeroDist= it->second;
		//Too much REprintf("Duplicate dist %g\n", xx);
	} else { //Build new node only when k==1
		ZeroDist=new CmpCluster();
		kendall[0.0]=ZeroDist;
	}
	TMRLOG(2,"Pearson loop and kendall tree\n");
	#ifdef DISTOBJECT_VERSION
	ij =-nelements;
	#endif
	for(i=0; i <nelements; i++) {
		iclustIndex = clusterid[i];
		sizes[iclustIndex]+=weights[i];
		#ifdef DISTMATRIX_VERSION
			ij=i*nelements;
		#else 
			// i_indiv=i+1;
			//base_indice=n*(i_indiv-1) - i_indiv*(i_indiv-1)/2 -i_indiv-1;
			// To be optimized
			ij += nelements-i-1;
			//ij  = i*nelements-((i+1)*i)/2-i-1;
		#endif
		
		if(weights[i]>0){
			// Take to not forget the diagonal!!!!
			ww = weights[i]*weights[i];
			// i obs is located in the same cluster as i :)
			wy+=ww;
			//Distance is 0 by definition
			ZeroDist->clustDist0+=ww;
			totweights+=ww;
			for(j=i+1;j <nelements;j++){
				if(weights[j]>0){
					ww = 2*weights[i]*weights[j];
					xx = distmatrix[ij+j];
					it=kendall.find(xx);
					if (it!=kendall.end()) {
						cmpclust= it->second;
						//Too much REprintf("Duplicate dist %g\n", xx);
					} else { //Build new node only when k==1
						cmpclust=new CmpCluster();
						kendall[xx]=cmpclust;
					}
					xw=ww*xx;
					xxw=xw*xx;
					wx+=xw;
					wx2+=xxw;
					if(clusterid[i] == clusterid[j]){
						errors[iclustIndex]+=xw;
						errors2[iclustIndex]+=xxw;
						wxy+=xw;
						//wxy2+=xxw;
						wy+=ww;
						cmpclust->clustDist0+=ww;
					}else{
						cmpclust->clustDist1+=ww;
					}
					
					totweights+=ww;
				}
			}
		}
	}
	xb = wx/totweights;
	yb = wy/totweights;
	covx = wx2/totweights - xb*xb;
	covy = wy/totweights - yb*yb;
	covxy = wxy/totweights - yb*xb;
	pearson = covxy/(R_pow(covx*covy, 0.5));
	stats[ClusterQualHPG] = -1.0*pearson; //HPG
	double nc=0, nd=0, currentclustdist0=0, currentclustdist1=0;
	double totdist0=wy, totdist1=totweights-wy, ntiesdist=0;
	double Smin=0, wSmin=wy, Smax=0, wSmax=totdist1, currentww=0;
	TMRLOG(2,"Kendall compute on %d values\n", kendall.size());
	for (it = kendall.begin();it != kendall.end();it++) {
        cmpclust=it->second;
		ww=cmpclust->clustDist1+cmpclust->clustDist0;
		if(ww>0) {
			if(currentww<=wSmin){
				if(currentww+ww>wSmin){
					Smin+=(wSmin-currentww)* it->first;
				}else{
					Smin+=ww* it->first;
				}
			}
			currentww+=ww;
			if(currentww>wSmax){
				if(currentww-ww<wSmax){
					Smax+=(currentww-wSmax)* it->first;
				}else{
					Smax+=ww*it->first;
				}
			}
			ntiesdist+=cmpclust->clustDist1*cmpclust->clustDist0;
			//Bottom of table
			nc+=cmpclust->clustDist1*currentclustdist0; //We have one and lesser distance have zero
			nd+=cmpclust->clustDist0*currentclustdist1; //We have zero and lesser distance have one
			//Adjusting
			currentclustdist0+=cmpclust->clustDist0;
			currentclustdist1+=cmpclust->clustDist1;
			
			//Up of table
			nc+= cmpclust->clustDist0*(totdist1-currentclustdist1);//We have zero and higher distance have one
			nd+= cmpclust->clustDist1*(totdist0-currentclustdist0); //We have one and higher distance have zero
		}
		
    }
	//Gamma ignoring ties
	stats[ClusterQualHG] =(nc-nd)/(nc+nd);
	//Somers D
	stats[ClusterQualHGSD] =(nc-nd)/(nc+nd+ntiesdist);
	stats[ClusterQualHC] =(wxy-Smin)/(Smax-Smin);
	double SSres=0;
	double SS2res=0;
	totweights=0;
	TMRLOG(2,"SS computation\n");
	for(i=0;i <nclusters;i++){
		//REprintf("TBS %d, sizes %f, errors %f errors2 %f\n", tablesizes[i], sizes[i], errors[i], errors2[i]);
		SSres+=errors[i]/sizes[i];
		SS2res+=errors2[i]/sizes[i];
		totweights+=sizes[i];
	}
	
	double SSexpl = wx/totweights-SSres;
	double SS2expl = wx2/totweights-SS2res;
	//REprintf("SSres %f, SS2res %f, SSexpl %f, SS2expl %f, Tot %f\n", SSres, SS2res, SSexpl, SS2expl, wx/totweights);
	double dncluster= (double)nclusters;
	stats[ClusterQualF] = (SSexpl/(dncluster-1.0))/(SSres/(totweights-dncluster)); //F
	stats[ClusterQualR] = (SSexpl/(SSres+SSexpl)); //R
	stats[ClusterQualF2] = (SS2expl/(dncluster-1.0))/(SS2res/(totweights-dncluster)); //F2
	
	stats[ClusterQualR2] = (SS2expl/(SS2res+SS2expl)); //R2
	
	//Computing ASW
	TMRLOG(2,"ASW statitstics\n");
	double asw_i=0;
	double asw_w=0;
	for (j=0;j<nclusters;j++) {
		errors2[j]=0.0;
		errors2[j+nclusters]=0.0;
	}
	for(i=0; i <nelements; i++) {
		if(weights[i]>0){
			iclustIndex = clusterid[i];
			double aik =0;
			double aik_w =0;
			double sik_i;
			double sik_w;
			for (j=0;j<nclusters;j++) {
				errors[j]=0.0;
			}
			#ifdef DISTMATRIX_VERSION
				ij=i*nelements;
				for(j=0;j <nelements;j++){
					if (i==j) {
						continue;
					}
					if(iclustIndex == clusterid[j]) {
						aik += weights[j]*distmatrix[ij+j];
					}
					else{
						errors[clusterid[j]]+=weights[j]*distmatrix[ij+j];
					}
				}
			#else 
				// i_indiv=i+1;
				//base_indice=n*(i_indiv-1) - i_indiv*(i_indiv-1)/2 -i_indiv-1;
				// To be optimized
				ij =DL_FIRST_INIT(i, nelements);
				for(j=0;j <i;j++){
					ij += DL_FIRST_INC(i, j, nelements);
					if(iclustIndex == clusterid[j]) {
						aik += weights[j]*distmatrix[DL_FIRST_ACCESS(ij, i, j)];
					}
					else{
						errors[clusterid[j]]+=weights[j]*distmatrix[DL_FIRST_ACCESS(ij, i, j)];
					}
				}
				ij = DL_SEC_INIT(i, nelements);
				for(j=i+1;j <nelements;j++){
					if(iclustIndex == clusterid[j]) {
						aik += weights[j]*distmatrix[DL_SEC_ACCESS(ij, i, j)];
					}
					else{
						errors[clusterid[j]]+=weights[j]*distmatrix[DL_SEC_ACCESS(ij, i, j)];
					}
				}
			#endif
			
			double bik =DBL_MAX;
			for (j=0; j<nclusters; j++) {
				if(j!=iclustIndex){
					if(bik>=(errors[j]/sizes[j])){
						bik=(errors[j]/sizes[j]);
					}
				}
			}
			//Avoid division by zero if  (sizes[iclustIndex]==weights[i]) (one observation per cluster)
			aik_w = aik/(sizes[iclustIndex]); 
			if(sizes[iclustIndex]<=1.0){
				aik =0;
			} else {
				aik /= (sizes[iclustIndex]-1);
			}
			
			sik_i = weights[i] * ((bik - aik)/fmax2(aik, bik));
			sik_w = weights[i] * ((bik - aik_w)/fmax2(aik_w, bik));
			//REprintf("aik %f, bik %f, sik %f\n", aik, bik, sik/weights[i]);
			errors2[iclustIndex]+=sik_i;
			errors2[iclustIndex+nclusters]+=sik_w;
			asw_i+=sik_i;
			asw_w+=sik_w;
		}
	}
	for (j=0;j<nclusters;j++) {
		errors2[j]/=sizes[j];
		errors2[j+nclusters]/=sizes[j];
	}
	stats[ClusterQualASWi] = asw_i/totweights; //R2
	stats[ClusterQualASWw] = asw_w/totweights; //R2
	return;
}

void INDIV_ASW_FUNCNAME(double * distmatrix, int * clusterid, double *weights, int nelements, int nclusters, double * asw_i, double * asw_w){
	TMRLOG(2,"Computing statitstics\n");
	int i, j, ij, iclustIndex;
	double *othergroups = (double*) R_alloc(nclusters, sizeof(double));
	double *sizes = (double*) R_alloc(nclusters, sizeof(double));
	for(i=0; i<nclusters; i++){
		othergroups[i]=0;
		sizes[i]=0;
	}
	
	//Compute cluster weighted sizes
	TMRLOG(2,"Pearson loop and kendall tree\n");
	for(i=0; i < nelements; i++) {
		sizes[clusterid[i]]+=weights[i];
	}
	
	for(i=0; i <nelements; i++) {
		iclustIndex = clusterid[i];
		double aik =0;
		double aik_w =0;
		for (j=0;j<nclusters;j++) {
			othergroups[j]=0.0;
		}
		
		#ifdef DISTMATRIX_VERSION
			ij=i*nelements;
			for(j=0;j <nelements;j++){
				if (i==j) {
					continue;
				}
				if(iclustIndex == clusterid[j]) {
					aik += weights[j]*distmatrix[ij+j];
				}
				else{
					othergroups[clusterid[j]]+=weights[j]*distmatrix[ij+j];
				}
			}
		#else 
			// i_indiv=i+1;
			//base_indice=n*(i_indiv-1) - i_indiv*(i_indiv-1)/2 -i_indiv-1;
			// To be optimized
			ij =DL_FIRST_INIT(i, nelements);
			for(j=0;j <i;j++){
				ij += DL_FIRST_INC(i, j, nelements);
				if(iclustIndex == clusterid[j]) {
					aik += weights[j]*distmatrix[DL_FIRST_ACCESS(ij, i, j)];
				}
				else{
					othergroups[clusterid[j]]+=weights[j]*distmatrix[DL_FIRST_ACCESS(ij, i, j)];
				}
			}
			ij = DL_SEC_INIT(i, nelements);
			for(j=i+1;j <nelements;j++){
				if(iclustIndex == clusterid[j]) {
					aik += weights[j]*distmatrix[DL_SEC_ACCESS(ij, i, j)];
				}
				else{
					othergroups[clusterid[j]]+=weights[j]*distmatrix[DL_SEC_ACCESS(ij, i, j)];
				}
			}
		#endif
		double bik =DBL_MAX;
		for (j=0; j<nclusters; j++) {
			if(j!=iclustIndex){
				if(bik>=(othergroups[j]/sizes[j])){
					bik=(othergroups[j]/sizes[j]);
				}
			}
		}
		aik_w = aik/ (sizes[iclustIndex]); 
		if(sizes[iclustIndex]<=1.0){
			aik =0;
		} else {
			aik /= (sizes[iclustIndex]-1);
		}
		
		asw_i[i] = ((bik - aik)/fmax2(aik, bik));
		asw_w[i] = ((bik - aik_w)/fmax2(aik_w, bik));
	}
}

void CLUSTERQUALITYSIMPLE_FUNCNAME(double * distmatrix, int * clusterid, double *weights, int nelements, double* stats, int nclusters, double * errors2){
	TMRLOG(2,"Computing statitstics\n");
	//double totweights=0, wxy=0,wxy2=0, wx=0, wy=0, wx2=0, ww, xx, covxy, covx, covy, pearson, xb, yb, xw, xxw;
	double totweights=0, wxy=0,wx=0, wy=0, wx2=0, ww, xx, covxy, covx, covy, pearson, xb, yb, xw, xxw;
	int i, j, ij, iclustIndex;
	double *errors = (double*) R_alloc(nclusters,sizeof(double));
	int *tablesizes = (int*) R_alloc(nclusters,sizeof(int));
	double *sizes = (double*) R_alloc(nclusters,sizeof(double));
	for(i=0;i <nclusters;i++){
		tablesizes[i]=-1;
		errors2[i]=0;
		errors[i]=0;
		sizes[i]=0;
	}
	
	//Maybe use  R_qsort_I to sort distances between individuals first
	//CmpCluster * cmpclust=NULL;
	TMRLOG(2,"Pearson loop and kendall tree\n");
	#ifdef DISTOBJECT_VERSION
		ij =-nelements;
	#endif
	for(i=0; i <nelements; i++) {
		iclustIndex = clusterid[i];
		sizes[iclustIndex]+=weights[i];
		#ifdef DISTMATRIX_VERSION
			ij=i*nelements;
		#else 
			// i_indiv=i+1;
			//base_indice=n*(i_indiv-1) - i_indiv*(i_indiv-1)/2 -i_indiv-1;
			// To be optimized
			ij += nelements-i-1;
		#endif
		
		if(weights[i]>0){
			// Take to not forget the diagonal!!!!
			ww = weights[i]*weights[i];
			// i obs is located in the same cluster as i :)
			wy+=ww;
			//Distance is 0 by definition
			totweights+=ww;

			for(j=i+1;j <nelements;j++){
				if(weights[j]>0){
					ww = weights[i]*weights[j];
					xx = distmatrix[ij+j];
					xw=ww*xx;
					xxw=xw*xx;
					wx+=xw;
					wx2+=xxw;
					if(clusterid[i] == clusterid[j]){
						errors[iclustIndex]+=xw;
						errors2[iclustIndex]+=xxw;
						wxy+=xw;
						//wxy2+=xxw;
						wy+=ww;
					}
					totweights+=ww;
				}
			}
		}
	}
	xb = wx/totweights;
	yb = wy/totweights;
	covx = wx2/totweights - xb*xb;
	covy = wy/totweights - yb*yb;
	covxy = wxy/totweights - yb*xb;
	pearson = covxy/(R_pow(covx*covy, 0.5));
	stats[ClusterQualHPG] = -1.0*pearson; //HPG
	double SSres=0;
	double SS2res=0;
	totweights=0;
	TMRLOG(2,"SS computation\n");
	for(i=0;i <nclusters;i++){
		//REprintf("TBS %d, sizes %f, errors %f errors2 %f\n", tablesizes[i], sizes[i], errors[i], errors2[i]);
		SSres+=errors[i]/sizes[i];
		SS2res+=errors2[i]/sizes[i];
		totweights+=sizes[i];
	}
	
	double SSexpl = wx/totweights-SSres;
	double SS2expl = wx2/totweights-SS2res;
	//REprintf("SSres %f, SS2res %f, SSexpl %f, SS2expl %f, Tot %f\n", SSres, SS2res, SSexpl, SS2expl, wx/totweights);
	double dncluster= (double)nclusters;
	stats[ClusterQualF] = (SSexpl/(dncluster-1.0))/(SSres/(totweights-dncluster)); //F
	stats[ClusterQualR] = (SSexpl/(SSres+SSexpl)); //R
	stats[ClusterQualF2] = (SS2expl/(dncluster-1.0))/(SS2res/(totweights-dncluster)); //F2
	stats[ClusterQualR2] = (SS2expl/(SS2res+SS2expl)); //R2
	
	return;
}

#endif //#ifdef CLUSTERQUALITY_INCLUDED
