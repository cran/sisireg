#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include <stdlib.h> // for NULL


extern void ssrC(int* funk, double* y, double* mu, int* anz, int* fn, int* ps, int* simanz);
extern void ssr_neC(double* x, double* y, double* mu, int* anz, int* fn, int* simanz);

static R_NativePrimitiveArgType ssr_t[] = {
  INTSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP
};
static R_NativePrimitiveArgType ssr_ne_t[] = {
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP
};


static const R_CMethodDef cMethods[] = {
  {"ssrC", (DL_FUNC) &ssrC, 7, ssr_t},
  {"ssr_neC", (DL_FUNC) &ssr_neC, 6, ssr_ne_t},
  {NULL, NULL, 0, NULL}
};

void R_init_sisireg(DllInfo* info) {
  R_registerRoutines(info, cMethods, NULL, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
  R_forceSymbols(info, TRUE);
}

