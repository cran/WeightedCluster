#include "PAM.h"
#include "PAMonce.h"
#include "kmedoid.h"
#include "kmedoidbase.h"
#include "cluster.h"
//	centroids = (int*) R_alloc(nclusters,sizeof(int));

int TRAMINER_DEBUG_LEVEL=0;
void finalizeKMedoidBase(SEXP ptr){
	KMedoidBase * sdo;
	sdo= static_cast<KMedoidBase *>(R_ExternalPtrAddr(ptr));
	delete sdo;
}


inline SEXP KMedoidBaseWorker(KMedoidBase *ds) {
    SEXP SDO, classname;
	PROTECT(classname = allocVector(STRSXP, 1));
	SET_STRING_ELT(classname, 0, mkChar("KMedoidBase"));
    SDO = R_MakeExternalPtr(ds, R_NilValue, R_NilValue);
    R_RegisterCFinalizerEx(SDO, (R_CFinalizer_t) finalizeKMedoidBase, TRUE);
    classgets(SDO, classname);
	UNPROTECT(1);
    return SDO;
}

extern "C" {
	SEXP RKmedoids(SEXP Snelement, SEXP diss, SEXP expr, SEXP rho, SEXP Scentroids, SEXP Snpass, SEXP Sweights, SEXP SclustMethod, SEXP SDebug, SEXP Sisdist){
		int old_debug=TRAMINER_DEBUG_LEVEL;
		TRAMINER_DEBUG_LEVEL=INTEGER(SDebug)[0];
		KMedoidBase * km;
		TMRLOG(1, "Starting clustering\n");
		if(INTEGER(SclustMethod)[0]==1){
			km = new KMedoid(Snelement, diss, expr,  rho, Scentroids, Snpass, Sweights, Sisdist);
		}else if(INTEGER(SclustMethod)[0]==2){
			km = new PAM(Snelement, diss, expr, rho, Scentroids, Snpass, Sweights, Sisdist);
		}else{
			km = new PAMonce(Snelement, diss, expr, rho, Scentroids, Snpass, Sweights, Sisdist);
		}
		SEXP ClusterAlgo;
		PROTECT(ClusterAlgo= KMedoidBaseWorker(km));
		km->findCluster();
		km->clean();
		UNPROTECT(1);
		TRAMINER_DEBUG_LEVEL=old_debug;
		return(km->getClustering());
	}
}
