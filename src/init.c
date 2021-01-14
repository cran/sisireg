#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdlib.h> // for NULL


extern void ssrC(int* funk, double* y, double* mu, int* anz, int* fn, int* simanz);
extern void ssr_neC(double* x, double* y, double* mu, int* anz, int* fn, int* simanz);

SEXP ssr3dC(SEXP koord, SEXP nb, SEXP nb1, SEXP z, SEXP fn, SEXP iteranz);
SEXP ssrndC(SEXP koord, SEXP nb_ms, SEXP nb_dst, SEXP nb_ps, SEXP y, SEXP fn, SEXP iteranz);

static R_NativePrimitiveArgType ssr_t[] = {
  INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};
static R_NativePrimitiveArgType ssr_ne_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};


static const R_CMethodDef cMethods[] = {
  {"ssrC", (DL_FUNC) &ssrC, 6, ssr_t},
  {"ssr_neC", (DL_FUNC) &ssr_neC, 6, ssr_ne_t},
  {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[] = {
  {"ssr3dC", (DL_FUNC) &ssr3dC, 6},
  {"ssrndC", (DL_FUNC) &ssrndC, 7},
  {NULL, NULL, 0}
};

void R_init_sisireg(DllInfo* info) {
  R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}

