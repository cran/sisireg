/**
 * Sign-Simplicity-Regression
 */
#include <stdio.h> // for output
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

//#define DEBUG // Log-Ausgaben

#define sign(x) ((x > 0) - (x < 0))
#define max(x,y) ((x<y)? y : x)
#define min(x,y) ((x>y)? y : x)


/**
 * QSOR-algorithm equidistant
 */
void ssrC(int* funk, double* y, double* mu, int* anz, int* fn, int *ps, int* simanz) {

#ifdef DEBUG
  fprintf(stdout, "ssr with parameters: funk=%d, anz=%d, fn=%d, ps=%d, simanz=%d\n", *funk, *anz, *fn, *ps, *simanz); fflush(stdout);
#endif

  int l = *funk;
  int n = *anz;

  int k = (*ps==1)? 3.3+1.44*log(n) : *fn; // 0.95%-quantile of max. run length
  int h = (*ps==1)? *fn : 0;

  if (l == 1) {
    for (int s = 0; s < *simanz; s++) {
      bool chng = false;
      for (int i = 1; i < n-1; i++) {
        double oldval = mu[i];
        int oldsig = sign(y[i] - mu[i]);
        mu[i] = (mu[i-1] + mu[i+1])/2;
        int newsig = sign(y[i] - mu[i]);
        if (*ps == 1) {
          if ((i > k) && (i < n-k)) {
            if (oldsig != newsig) {
              int sum = 0;
              for (int m = -k; m <= k; m++)
                sum += sign(y[i+m] - mu[i+m]);
              if (abs(sum) > h) mu[i] = oldval;
            }
          }
        } else {
          for (int j = max(0,i-k); j < min(n-k, i+k); j++) {
            if (oldsig != newsig) {
              int sum = 0;
              for (int m = 0; m <= k; m++)
                sum += sign(y[j+m] - mu[j+m]);
              if (abs(sum) > k) {
                mu[i] = oldval;
                break;
              }
            }
          }
        }
        if (mu[i] != oldval) chng = true;
      }
      if (!chng) {
        break;
      }
    }
  } else if (l == 2) {
    // QSOR-matrix
    double a2 = sqrt(2.);
    double a10 = sqrt(10.);
    double a12 = sqrt(12.);
    double q12 = 1. / (a2 * a10);
    double q13 = 1. / (a2 * a12);
    double q23 = 1. / (a10 * a12);
    double a012 = (-4.) * q12;
    double a023 = (-8.) * q23;
    double a013 = 2. * q13;
    double a024 = 2. * q23;
    double as2 = 2. / (12.);
    double as8 = (-8.) / (12.);
    double omega = 1.9;

    // standardize the data
    mu[0] = mu[0] * a2;
    y[0] = y[0] * a2;
    mu[1] = mu[1] * a10;
    y[1] = y[1] * a10;
    for (int i = 2; i < n-2; i++) {
      mu[i] = mu[i] * a12;
      y[i] = y[i] * a12;
    }
    mu[n-2] = mu[n-2] * a10;
    y[n-2] = y[n-2] * a10;
    mu[n-1] = mu[n-1] * a2;
    y[n-1] = y[n-1] * a2;

#ifdef DEBUG
    fprintf(stdout, "calculating:["); for (int i=0;i<50;i++) {fprintf(stdout, ".");} fprintf(stdout, "]");
    for (int i=0;i<51;i++) {fprintf(stdout, "\b");} fflush(stdout);
#endif
    // iterations
    for (int s = 0; s < *simanz; s++) {
      bool chng = false;
#ifdef DEBUG
      if (s % (max(*simanz,50) / 50) == 0) { fprintf(stdout, "X"); fflush(stdout);}
#endif
      for (int i = 1; i < n-1; i++) {
        double oldval = mu[i];
        int oldsig = sign(y[i] - mu[i]);

        if (i == 1)   mu[i] = mu[i] - omega*(               mu[i-1]*a012 + mu[i] + mu[i+1]*a023 + mu[i+2]*a024);
        if (i == 2)   mu[i] = mu[i] - omega*(mu[i-2]*a013 + mu[i-1]*a023 + mu[i] + mu[i+1]*as8  + mu[i+2]*as2 );
        if (i == 3)   mu[i] = mu[i] - omega*(mu[i-2]*a024 + mu[i-1]*as8  + mu[i] + mu[i+1]*as8  + mu[i+2]*as2 );
        if ((i > 3) && (i < n-4))
                      mu[i] = mu[i] - omega*(mu[i-2]*as2  + mu[i-1]*as8  + mu[i] + mu[i+1]*as8  + mu[i+2]*as2 );
        if (i == n-4) mu[i] = mu[i] - omega*(mu[i-2]*as2  + mu[i-1]*as8  + mu[i] + mu[i+1]*as8  + mu[i+2]*a024);
        if (i == n-3) mu[i] = mu[i] - omega*(mu[i-2]*as2  + mu[i-1]*as8  + mu[i] + mu[i+1]*a023 + mu[i+2]*a013);
        if (i == n-2) mu[i] = mu[i] - omega*(mu[i-2]*a024 + mu[i-1]*a023 + mu[i] + mu[i+1]*a012);

        int newsig = sign(y[i] - mu[i]);
        if (*ps == 1) {
          if ((i > k) && (i < n-k)) {
            if (oldsig != newsig) {
              int sum = 0;
              for (int m = -k; m <= k; m++)
                sum += sign(y[i+m] - mu[i+m]);
              if (abs(sum) > h) mu[i] = oldval;
            }
          }
        } else {
          for (int j = max(0,i-k); j <= min(n-k, i+k); j++) {
            if (oldsig != newsig) {
              int sum = 0;
              for (int m = 0; m <= k; m++)
                sum += sign(y[j+m] - mu[j+m]);
              if (abs(sum) > k) {
                mu[i] = oldval;
                break;
              }
            }
          }
        }
        if (mu[i] != oldval) chng = true;
      }
      if (!chng) {
        break;
      }
    }
#ifdef DEBUG
    fprintf(stdout, "\n");
    fflush(stdout);
#endif
    // re-transform the data
    mu[0] = mu[0] / a2;
    y[0] = y[0] / a2;
    mu[1] = mu[1] / a10;
    y[1] = y[1] / a10;
    for (int i = 2; i < n-2; i++) {
      mu[i] = mu[i] / a12;
      y[i] = y[i] / a12;
    }
    mu[n-2] = mu[n-2] / a10;
    y[n-2] = y[n-2] / a10;
    mu[n-1] = mu[n-1] / a2;
    y[n-1] = y[n-1] / a2;
  }

}


/**
 * QSOR-algorithm L1 non-equidistant
 * x is presented in sorted order
 */
void ssr_neC(double* x, double* y, double* mu, int* anz, int* fn, int* simanz) {

  // parameters
  int n = *anz;
  int k = 3.3+1.44*log(n); // 95%-quantile of max. run length
  int h = *fn;

  int sim = *simanz;
  for (int j = 0; j < sim; j++) {
    bool chng = false;
    for (int i = 1; i < n-1; i++) {
      // save old value
      double oldval = mu[i];
      int oldsig = sign(y[i] - mu[i]);
      // equal values are permitted
      if (x[i-1] != x[i+1]) {
        mu[i] = mu[i-1] + (x[i]-x[i-1])*(mu[i+1] - mu[i-1])/(x[i+1]- x[i-1]);
      } else {
        mu[i] = (mu[i-1] + mu[i+1])/2;
      }
      // check partial sums
      int newsig = sign(y[i] - mu[i]);
      if ((i > k) && (i < n-k)) {
        if (oldsig != newsig) {
          int sum = 0;
          for (int m = -k; m <= k; m++)
            sum += sign(y[i+m] - mu[i+m]);
          if (abs(sum) > h) mu[i] = oldval;
        }
      }
      if (mu[i] != oldval) chng = true;
    }
    if (!chng) {
      break;
    }
  }
}

