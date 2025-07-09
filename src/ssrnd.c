/**
 * Sign-Simplicity-Regression: n-dimensional fUnctions
 */
#include <R.h>
#include <Rinternals.h>

#include <stdio.h> // for output
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>


//#define DEBUG // Log-Ausgaben

#define sign(x) ((x > 0) - (x < 0))
#define max(x,y) (x>y?x:y)


/**
 * chi-function: checks, if the newly calculated value is valid
 */
bool chi_nd(int** nb_vec, int n, int* nb_anz, int ind, double* z, double* mu, double fn) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < nb_anz[i]; j++) {
      if (nb_vec[i][j] == ind) {
        // sum of the residual signs
        int sum = 0;
        for (int m = 0; m < nb_anz[i]; m++) {
          sum += sign(z[nb_vec[i][m]] - mu[nb_vec[i][m]]);
        }
        if (abs(sum) > fn) {
          return FALSE;
        }
        break;
      }
    }
  }
  return TRUE;
}

/**
 * calculates the exponential-distance-weighted mean
 */
double wmean_nd(int* nb_vec[], double* nb_dst[], int* nb_anz, double* p_mu, int i) {
  double val = 0;
  double w = 0;
  // weighted mean
  for (int m = 0; m < nb_anz[i]; m++) {
    double d = nb_dst[i][m];
    if (d != 0) { // omit the reference point
      double g = exp(-d);
      w += g;
      val += (p_mu[nb_vec[i][m]]*g);
    }
  }
  return val / w;

}

/**
 * extracts the coordinates
 */
void getCoordinates(SEXP koord, int n, int dim, double* pd_koord) {
  // internal input vector
  for (int j = 0; j < dim; j++) {
    SEXP koord_vec_x = coerceVector(VECTOR_ELT(koord, j), REALSXP);
    double* pd_vec = REAL(koord_vec_x);
    // new internal output vector
    for (int i = 0; i < n; i++) pd_koord[j+i*dim] = pd_vec[i];
  }
}

/**
 * extracts the indices and the distances for the neighborhood
 */
void getNeighborhood(SEXP nb_ind, SEXP nb_dst, int n, int* pi_nb_vec[], int* i_nb_anz, double* pd_dst_vec[]) {
  for (int i = 0; i < n; i++) {
    // indices of neighborhood - internal input vector
    SEXP nbvec = coerceVector(VECTOR_ELT(nb_ind, i), INTSXP);
    int* pi_ind = INTEGER(nbvec);
    R_xlen_t m = xlength(nbvec);
    i_nb_anz[i] = m;
    // new internal output vector
    SEXP neighbours = PROTECT(allocVector(INTSXP, m));
    int* pi_neighbours = INTEGER(neighbours);
    for (R_xlen_t j = 0; j < m; j++) {
      pi_neighbours[j] = pi_ind[j] - 1; // C is zero-indexed
    }
    pi_nb_vec[i] = pi_neighbours;

    // distances of neighborhood
    if (nb_dst != NULL) {
      // internal input vector
      SEXP dstvec = coerceVector(VECTOR_ELT(nb_dst, i), REALSXP);
      double* pd_dst = REAL(dstvec);
      // new internal output vector
      SEXP dst_neighbours = PROTECT(allocVector(REALSXP, m));
      double* pd_dst_neighbours = REAL(dst_neighbours);
      for (R_xlen_t j = 0; j < m; j++) {
        pd_dst_neighbours[j] = pd_dst[j];
      }
      pd_dst_vec[i] = pd_dst_neighbours;
    }
  }
}

/**
 * calculates the minimal surface
 */
SEXP ssrndC(SEXP koord, SEXP nb_ms, SEXP nb_dst, SEXP nb_ps, SEXP y, SEXP fn, SEXP iteranz) {

  double d_fn = asReal(fn);
  int i_iteranz = asInteger(iteranz);

  // dimensions
  SEXP cols = coerceVector(koord, VECSXP);
  SEXP rows = coerceVector(VECTOR_ELT(cols, 1), REALSXP);
  int n = length(rows);
  int dim = length(cols);
#ifdef DEBUG
  fprintf(stdout, "ssrndC: %d observations in %d dimensions\n", n, dim); fflush(stdout);
#endif
  // coordinates
  double koord_vec[n][dim];
  // internal input vector
  getCoordinates(koord, n, dim, &(koord_vec[0][0]));

  // neighborhood for minimal surface
  int      i_nbms_anz[n];
  int*    pi_nbms_vec[n];
  double* pd_nbdst_vec[n];
  getNeighborhood(nb_ms, nb_dst, n, pi_nbms_vec, i_nbms_anz, pd_nbdst_vec);

  // neighborhood PS
  int   i_nbps_anz[n];
  int* pi_nbps_vec[n];
  getNeighborhood(nb_ps, NULL, n, pi_nbps_vec, i_nbps_anz, NULL);

  // create target vector
  double* p_y =  REAL(y);
  SEXP mu = PROTECT(allocVector(REALSXP, n));
  double* p_mu = REAL(mu);

  // start values: running mean of the neighborhood
  for (int i = 0; i < n; i++)  {
    p_mu[i] = wmean_nd(pi_nbms_vec, pd_nbdst_vec, i_nbms_anz, p_y, i);
  }
  for (int i = 0; i < n; i++) {
    if (!chi_nd(pi_nbps_vec, n, i_nbps_anz, i, p_y, p_mu, d_fn))  p_mu[i] = p_y[i];
  }


  bool change = TRUE;
  // iterations
#ifdef DEBUG
  fprintf(stdout, "ssrndC: calculating:["); for (int i=0;i<50;i++) {fprintf(stdout, ".");} fprintf(stdout, "]");
  for (int i=0;i<51;i++) {fprintf(stdout, "\b");} fflush(stdout);
#endif
  for (int s = 1; s <= i_iteranz; s++) {
    change = FALSE;
#ifdef DEBUG
    if (s % (max(i_iteranz,50) / 50) == 0) { fprintf(stdout, "X"); fflush(stdout);}
#endif
    for (int i = 0; i < n; i++) {
      double oldval = p_mu[i];
      double newval = wmean_nd(pi_nbms_vec, pd_nbdst_vec, i_nbms_anz, p_mu, i);
      p_mu[i] = newval;
      //if (!chi_nd(pi_nbps_vec, n, i_nbps_anz, i, p_y, p_mu, d_fn)) {
      //  if (sign(p_y[i] - newval) == sign(p_y[i] - oldval))   p_mu[i] = newval;
      //  else                                                  p_mu[i] = p_y[i];
      //}
      if (!chi_nd(pi_nbps_vec, n, i_nbps_anz, i, p_y, p_mu, d_fn)) {
        p_mu[i] = oldval; // intial mu fulfills ps criterion
      }
      if (p_mu[i] != oldval)  change = TRUE;
    }
    if (!change)  {
      break;
    }
  }
#ifdef DEBUG
  fprintf(stdout, "\n");
  fflush(stdout);
#endif
  UNPROTECT(3*n + 1);
  return mu;
}

