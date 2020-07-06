#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP RClusterQualBootSeveral(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RClusterComputeIndivASW(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RClusterQual(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RClusterQualKendall(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP RClusterQualKendallFactory();
extern SEXP RKmedoids(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"RClusterQualBootSeveral",    (DL_FUNC) &RClusterQualBootSeveral,    11},
  {"RClusterComputeIndivASW",    (DL_FUNC) &RClusterComputeIndivASW,     5},
  {"RClusterQual",               (DL_FUNC) &RClusterQual,                5},
  {"RClusterQualKendall",        (DL_FUNC) &RClusterQualKendall,         6},
  {"RClusterQualKendallFactory", (DL_FUNC) &RClusterQualKendallFactory,  0},
  {"RKmedoids",                  (DL_FUNC) &RKmedoids,                  10},
  {NULL, NULL, 0}
};

void R_init_WeightedCluster(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
