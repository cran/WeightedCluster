#include "cluster.h"
/*R_INLINE int getSizeIndex(int* tablesizes, int clust, int maxsize){
	int i=0;
	while(i<maxsize){
		if(tablesizes[i]==clust){
			return i;
		} else if(tablesizes[i]<0){
			tablesizes[i] = clust;
			return i;
		}
		i++;
	}
	error("Cluster not found in table\n");
	return i;
	
}*/

#include <map>


class CmpCluster{
	public:
	double clustDist0;
	double clustDist1;
	CmpCluster():clustDist0(0), clustDist1(0){}
	~CmpCluster(){}
	
};
#define ClusterQualHPG 0
#define ClusterQualHG 1
#define ClusterQualHGSD 2
#define ClusterQualASWi 3
#define ClusterQualASWw 4
#define ClusterQualF 5
#define ClusterQualR 6
#define ClusterQualF2 7
#define ClusterQualR2 8
#define ClusterQualHC 9
#define ClusterQualNumStat 10


typedef std::map<double, CmpCluster *> KendallTree;
typedef std::map<double, CmpCluster *>::iterator KendallTreeIterator;

void finalizeKendall(SEXP ptr){
	KendallTree * kendall;
	kendall= static_cast<KendallTree *>(R_ExternalPtrAddr(ptr));
	KendallTreeIterator it;
	for (it = kendall->begin();it != kendall->end();it++) {
		delete it->second;
	}
	delete kendall;
}

SEXP kendallFactory(KendallTree *kendall) {
    SEXP SDO, classname;
	PROTECT(classname = allocVector(STRSXP, 1));
	SET_STRING_ELT(classname, 0, mkChar("KendallTree"));
    SDO = R_MakeExternalPtr(kendall, R_NilValue, R_NilValue);
    R_RegisterCFinalizerEx(SDO, (R_CFinalizer_t) finalizeKendall, TRUE);
    classgets(SDO, classname);
	UNPROTECT(1);
    return SDO;
}


void resetKendallTree(KendallTree * kendall){
	TMRLOG(2, "Resetting kendall\n");
	KendallTreeIterator it;
	for (it = kendall->begin();it != kendall->end();it++) {
		it->second->clustDist0=0;
		it->second->clustDist1=0;
	}
}


#define CLUSTERQUALITY_INCLUDED
	//Including version based on dist object
	#define DISTOBJECT_VERSION
	#define CLUSTERQUALITY_FUNCNAME clusterquality_dist
	#define INDIV_ASW_FUNCNAME indiv_asw_dist
	#define CLUSTERQUALITYSIMPLE_FUNCNAME clusterqualitySimple_dist
		#include "clusterqualitybody.cpp"
	#undef DISTOBJECT_VERSION
	#undef CLUSTERQUALITY_FUNCNAME 
	#undef INDIV_ASW_FUNCNAME 
	#undef CLUSTERQUALITYSIMPLE_FUNCNAME

	//Including version based on a distance matrix
	#define DISTMATRIX_VERSION
	#define CLUSTERQUALITY_FUNCNAME clusterquality
	#define INDIV_ASW_FUNCNAME indiv_asw
	#define CLUSTERQUALITYSIMPLE_FUNCNAME clusterqualitySimple
		#include "clusterqualitybody.cpp"
	#undef DISTMATRIX_VERSION
	#undef CLUSTERQUALITY_FUNCNAME 
	#undef INDIV_ASW_FUNCNAME 
	#undef CLUSTERQUALITYSIMPLE_FUNCNAME

#undef CLUSTERQUALITY_INCLUDED

extern "C" {
	SEXP RClusterQual(SEXP diss, SEXP cluster, SEXP weightSS, SEXP numclust, SEXP isdist){
		int nclusters=INTEGER(numclust)[0];
		SEXP ans, stats, asw;
		PROTECT(ans = allocVector(VECSXP, 2));
		PROTECT(stats = allocVector(REALSXP, ClusterQualNumStat));
		PROTECT(asw = allocVector(REALSXP, 2*nclusters));
		SET_VECTOR_ELT(ans, 0, stats);
		SET_VECTOR_ELT(ans, 1, asw);
		KendallTree kendall;
		if(INTEGER(isdist)[0]){
			clusterquality_dist(REAL(diss), INTEGER(cluster), REAL(weightSS), length(cluster), REAL(stats), nclusters, REAL(asw), kendall);
		} else {
			clusterquality(REAL(diss), INTEGER(cluster), REAL(weightSS), length(cluster), REAL(stats), nclusters, REAL(asw), kendall);
		}
		KendallTreeIterator it;
		for (it = kendall.begin();it != kendall.end();it++) {
			delete it->second;
		}
		UNPROTECT(3);
		return ans;
		
	}
	
	SEXP RClusterComputeIndivASW(SEXP diss, SEXP cluster, SEXP weightSS, SEXP numclust, SEXP isdist){
		int nclusters=asInteger(numclust);
		SEXP ans, asw_i, asw_w;
		PROTECT(asw_i = allocVector(REALSXP, length(cluster)));
		PROTECT(asw_w = allocVector(REALSXP, length(cluster)));
		PROTECT(ans = allocVector(VECSXP, 2));
		SET_VECTOR_ELT(ans, 0, asw_i);
		SET_VECTOR_ELT(ans, 1, asw_w);
		if(INTEGER(isdist)[0]){
			indiv_asw_dist(REAL(diss), INTEGER(cluster), REAL(weightSS), length(cluster), nclusters, REAL(asw_i), REAL(asw_w));
		} else {
			indiv_asw(REAL(diss), INTEGER(cluster), REAL(weightSS), length(cluster), nclusters, REAL(asw_i), REAL(asw_w));
		}
		UNPROTECT(3);
		return ans;
		
	}
	SEXP RClusterQualSimple(SEXP diss, SEXP cluster, SEXP weightSS, SEXP numclust, SEXP isdist){
		int nclusters=INTEGER(numclust)[0];
		SEXP stats, asw;
		PROTECT(stats = allocVector(REALSXP, ClusterQualNumStat));
		PROTECT(asw = allocVector(REALSXP, nclusters));
		if(INTEGER(isdist)[0]){
			clusterqualitySimple_dist(REAL(diss), INTEGER(cluster), REAL(weightSS), length(cluster), REAL(stats), nclusters, REAL(asw));
		} else {
			clusterqualitySimple(REAL(diss), INTEGER(cluster), REAL(weightSS), length(cluster), REAL(stats), nclusters, REAL(asw));
		}
		UNPROTECT(2);
		return stats;
		
	}
	
	
	SEXP RClusterQualSimpleBoot(SEXP diss, SEXP cluster, SEXP weightSS, SEXP numclust, SEXP Rs, SEXP expr, SEXP rho, SEXP nindivS, SEXP isdist){
		int nclusters=INTEGER(numclust)[0];
		int ncase = length(cluster);
		int R = INTEGER(Rs)[0];
		int nindiv = INTEGER(nindivS)[0];
		int simplestat[5] = {ClusterQualHPG, ClusterQualF, ClusterQualR, ClusterQualF2, ClusterQualR2};
		int i, r;
		double * weights = new double[ncase];
		
		double * stat=new double[ClusterQualNumStat];
		double *asw= new double[nclusters];
		SEXP randomSample, st;
		int * rs;
		PROTECT(st = allocMatrix(REALSXP, R, 5));
		double *stt =REAL(st);
		if(INTEGER(isdist)[0]){
			clusterqualitySimple_dist(REAL(diss), INTEGER(cluster), REAL(weightSS), length(cluster), stat, nclusters, asw);
		} else {
			clusterqualitySimple(REAL(diss), INTEGER(cluster), REAL(weightSS), length(cluster), stat, nclusters, asw);
		}
		for(i=0; i<5;i++){
			stt[i*R] = stat[simplestat[i]];
		}
		for(r=1; r<R; r++){
			PROTECT(randomSample = eval(expr, rho));
			rs=INTEGER(randomSample);
			for(i=0; i<ncase; i++){	
				weights[i]=0;
			}
			for(i=0; i<nindiv; i++){	
				weights[rs[i]]++;
			}
			UNPROTECT(1);
			if(INTEGER(isdist)[0]){
				clusterqualitySimple_dist(REAL(diss), INTEGER(cluster), weights, length(cluster), stat, nclusters, asw);
			} else {
				clusterqualitySimple(REAL(diss), INTEGER(cluster), weights, length(cluster), stat, nclusters, asw);
			}
			for(i=0; i<5;i++){
				stt[r+i*R] = stat[simplestat[i]];
			}
			
		}
		//ww <- tabulate(sample.int(nrow(diss), size=totweights, replace=TRUE, prob=prob), nbins=nrow(diss))
		delete [] stat;
		delete [] asw;
		UNPROTECT(1);
		return st;
		
	}
	SEXP RClusterQualInitBoot(){
		return(kendallFactory(new KendallTree()));
	}
	SEXP RClusterQualBoot(SEXP diss, SEXP cluster, SEXP weightSS, SEXP numclust, SEXP kendallS, SEXP isdist){
		int nclusters=INTEGER(numclust)[0];
		SEXP ans, stats, asw;
		PROTECT(ans = allocVector(VECSXP, 2));
		PROTECT(stats = allocVector(REALSXP, ClusterQualNumStat));
		PROTECT(asw = allocVector(REALSXP, 2*nclusters));
		SET_VECTOR_ELT(ans, 0, stats);
		SET_VECTOR_ELT(ans, 1, asw);
		KendallTree * kendall;
		kendall= static_cast<KendallTree *>(R_ExternalPtrAddr(kendallS));
		resetKendallTree(kendall);
		if(INTEGER(isdist)[0]){
			clusterquality_dist(REAL(diss), INTEGER(cluster), REAL(weightSS), length(cluster), REAL(stats), nclusters, REAL(asw), (*kendall));
		} else {
			clusterquality(REAL(diss), INTEGER(cluster), REAL(weightSS), length(cluster), REAL(stats), nclusters, REAL(asw), (*kendall));
		}
		UNPROTECT(3);
		return ans;
		
	}
}
