/**
 * Sign-Simplicity-Regression: 3D fUnctions
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

// euclidian distance in two dimensions
#define dist(x,y) (sqrt((x[0]-y[0])*(x[0]-y[0]) + (x[1]-y[1])*(x[1]-y[1])))

/**
 * chi-function: checks, if the newly calculated value is valid
 */
bool chi(int** nb_vec, int n, int* nb_anz, int ind, double* z, double* mu, double fn) {
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
 * calculates the reciprocal-distance-weighted mean
 */
double wmean(double koord_vec[][2], int* nb1_vec[], int* nb1_anz, double* p_mu, int i) {
  double val = 0;
  double w = 0;
  // weighted mean
  for (int m = 0; m < nb1_anz[i]; m++) {
    double d = dist(koord_vec[i], koord_vec[nb1_vec[i][m]]);
    if (d != 0) { // omit the reference point
      w += 1/d;
      val += (p_mu[nb1_vec[i][m]]/d);
    }
  }
  return val / w;

}


/**
 * calculates the exponential-distance-weighted mean
 */
double wmean_exp(double koord_vec[][2], int* nb1_vec[], int* nb1_anz, double* p_mu, int i) {
  double val = 0;
  double w = 0;
  // weighted mean
  for (int m = 0; m < nb1_anz[i]; m++) {
    double d = dist(koord_vec[i], koord_vec[nb1_vec[i][m]]);
    if (d != 0) { // omit the reference point
      double g = exp(-d);
      w += g;
      val += (p_mu[nb1_vec[i][m]]*g);
    }
  }
  return val / w;

}


/**
 * calculates the minimal surface
 */
SEXP ssr3dC(SEXP koord, SEXP nb, SEXP nb1, SEXP z, SEXP fn, SEXP iteranz) {

  nb = coerceVector(nb, VECSXP);
  int n = length(nb);
  double d_fn = asReal(fn);
  int i_iteranz = asInteger(iteranz);

  // coordinates
  double koord_vec[n][2];
  // internal input vector
  for (int j = 0; j < 2; j++) {
    SEXP koord_vec_x = coerceVector(VECTOR_ELT(koord, j), REALSXP);
    double* pd_vec = REAL(koord_vec_x);
    // new internal output vector
    for (int i = 0; i < n; i++) {
      koord_vec[i][j] = pd_vec[i];
    }
  }

  // neighborhood PS
  int nb_anz[n];
  int* nb_vec[n];
  for (int i = 0; i < n; i++) {
    // internal input vector
    SEXP nbvec = coerceVector(VECTOR_ELT(nb, i), INTSXP);
    int* p_vec = INTEGER(nbvec);
    R_xlen_t m = xlength(nbvec);
    nb_anz[i] = m;
    // new internal output vector
    SEXP neighbours = PROTECT(allocVector(INTSXP, m));
    int* p_neighbours = INTEGER(neighbours);
    for (R_xlen_t j = 0; j < m; j++) {
      p_neighbours[j] = p_vec[j] - 1; // C is zero-indexed
    }
    nb_vec[i] = p_neighbours;
  }

  // neighborhood for mean value
  int nb1_anz[n];
  int* nb1_vec[n];
  for (int i = 0; i < n; i++) {
    // internal input vector
    SEXP nb1vec = coerceVector(VECTOR_ELT(nb1, i), INTSXP);
    int* p_vec1 = INTEGER(nb1vec);
    R_xlen_t m1 = xlength(nb1vec);
    nb1_anz[i] = m1;
    // new internal output vector
    SEXP neighbours1 = PROTECT(allocVector(INTSXP, m1));
    int* p_neighbours1 = INTEGER(neighbours1);
    for (R_xlen_t j = 0; j < m1; j++) {
      p_neighbours1[j] = p_vec1[j] - 1; // C is zero-indexed
    }
    nb1_vec[i] = p_neighbours1;
  }

  // create target vector
  double* p_z =  REAL(z);
  SEXP mu = PROTECT(allocVector(REALSXP, n));
  double* p_mu = REAL(mu);


  // start values: running mean of the neighborhood
  for (int i = 0; i < n; i++)  {
      p_mu[i] = wmean_exp(koord_vec, nb_vec, nb_anz, p_z, i);
  }
  for (int i = 0; i < n; i++) {
    if (!chi(nb_vec, n, nb_anz, i, p_z, p_mu, d_fn))  p_mu[i] = p_z[i];
  }


  bool change = TRUE;
  // iterations
#ifdef DEBUG
  fprintf(stdout, "calculating:["); for (int i=0;i<50;i++) {fprintf(stdout, ".");} fprintf(stdout, "]");
  for (int i=0;i<51;i++) {fprintf(stdout, "\b");} fflush(stdout);
#endif
  for (int s = 1; s <= i_iteranz; s++) {
    change = FALSE;
#ifdef DEBUG
    if (s % (max(i_iteranz,50) / 50) == 0) { fprintf(stdout, "X"); fflush(stdout);}
#endif
    for (int i = 0; i < n; i++) {
      double oldval = p_mu[i];
      double newval = wmean_exp(koord_vec, nb1_vec, nb1_anz, p_mu, i);
      p_mu[i] = newval;
      if (!chi(nb_vec, n, nb_anz, i, p_z, p_mu, d_fn)) {
        if (sign(p_z[i] - newval) == sign(p_z[i] - oldval))   p_mu[i] = newval;
        else                                                  p_mu[i] = p_z[i];
      }
      if (fabs(p_mu[i]-oldval) > 0.0000001)  change = TRUE;
    }
    if (!change)  {
      break;
    }
  }
#ifdef DEBUG
  fprintf(stdout, "\n");
  fflush(stdout);
#endif
  UNPROTECT(2*n + 1);
  return mu;
}


