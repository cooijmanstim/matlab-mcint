#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>

#include "mex.h"
#include "matrix.h"

/* for windows */
#ifndef snprintf
#define snprintf _snprintf
#endif

char *getNewString(const mxArray *ms) {
  size_t len;
  char *s;
  
  len = mxGetNumberOfElements(ms);
  s = malloc(len * sizeof(char) + 1);
  mxGetString(ms, s, len + 1);
  return s;
}

/* If an exception occurs in the MATLAB function, MATLAB's default
 * behavior seems to be to crash.  Let's try to deal with errors somewhat
 * properly.
 * This is also used for GSL errors.
 */
mxArray *exception = NULL;

void gslErrorHandler(const char *reason, const char *file, int line, int gsl_errno) {
  char *msgIdent, *msgString;
  const char *msgIdentPrefix, *msgStringPrefix;
  size_t msgIdentLength, msgStringLength;
  mxArray *rhs[2];
  mxArray *ex1, *ex2;

  // hate hate hate
  msgIdentPrefix = "MCI:GSL:errno";
  msgIdentLength = 20 + strlen(msgIdentPrefix);
  msgIdent = malloc(msgIdentLength * sizeof(char));
  snprintf(msgIdent, msgIdentLength, "%s%i", msgIdentPrefix, gsl_errno);
  rhs[0] = mxCreateString(msgIdent);

  msgStringPrefix = reason;
  msgStringLength = 20 + strlen(file) + strlen(msgStringPrefix);
  msgString = malloc(msgStringLength * sizeof(char));
  snprintf(msgString, msgStringLength, "'%s' at %s:%i", msgStringPrefix, file, line);
  rhs[1] = mxCreateString(msgString);

  ex2 = mexCallMATLABWithTrap(1, &ex1, 2, rhs, "MException");
  exception = ex2 == NULL ? ex1 : ex2;
}

double integrand(double x[], size_t dim, void *p) {
  mxArray *lhs[1], *rhs[2];
  double *x2;
  int i;
  double y;
  
  /* If an exception has occurred, don't call the MATLAB function again.
   * If there were a way to stop the GSL routines' continuing calls to the
   * present function, here would be a good place to do so. */
  if (exception != NULL)
      return 0;
  
  rhs[0] = *((mxArray **)p);
  rhs[1] = mxCreateDoubleMatrix(dim, 1, mxREAL);

  assert(mxIsClass(rhs[0], "function_handle"));

  /* We could do this without copying, but I'm not sure how safe it is to
   * change a Matlab object's data pointer to point to a region of memory
   * that's owned by the GSL routines.  It might be freed or the data might
   * be changed by GSL without Matlab knowing it. */
  x2 = mxGetPr(rhs[1]);
  for (i = 0; i < dim; i++)
    x2[i] = x[i];

  exception = mexCallMATLABWithTrap(1, lhs, 2, rhs, "feval");
  
  /* lhs is probably corrupt here, so don't use it. */
  if (exception != NULL)
      return 0;
  
  y = mxGetScalar(lhs[0]);

  mxDestroyArray(rhs[1]);
  mxDestroyArray(lhs[0]);
  return y;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const int minnrhs = 6;
  char *algorithm;
  int dim;
  double *A, *B;
  const mxArray *mxf;
  int calls;
  gsl_monte_function gslf;
  gsl_rng *rng;
  double result = 0, abserr = 0;
  int j;

  gsl_set_error_handler(&gslErrorHandler);
  
  if (nrhs < minnrhs)
    mexErrMsgIdAndTxt("MCI:BadArgument", "not enough arguments given");

  if (!mxIsChar(prhs[0]))
    mexErrMsgIdAndTxt("MCI:BadArgument", "algorithm must be a string; use one of {plain,miser,vegas}");
  algorithm = getNewString(prhs[0]);

  if (!mxIsDouble(prhs[1]) || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1)
    mexErrMsgIdAndTxt("MCI:BadArgument", "dim must be a nonnegative integer");
  dim = (int)mxGetScalar(prhs[1]);

  if (!mxIsDouble(prhs[2]) || mxGetNumberOfElements(prhs[2]) != dim)
    mexErrMsgIdAndTxt("MCI:BadArgument", "A must be a vector of length dim specifying the lower bounds of each component of x");
  if (!mxIsDouble(prhs[3]) || mxGetNumberOfElements(prhs[3]) != dim)
    mexErrMsgIdAndTxt("MCI:BadArgument", "B must be a vector of length dim specifying the upper bounds of each component of x");

  A = malloc(sizeof(double) * dim);
  for (j = 0; j < dim; j++)
    A[j] = mxGetPr(prhs[2])[j];

  B = malloc(sizeof(double) * dim);
  for (j = 0; j < dim; j++)
    B[j] = mxGetPr(prhs[3])[j];

  if (!mxIsClass(prhs[4], "function_handle"))
    mexErrMsgIdAndTxt("MCI:BadArgument", "f must be a function handle");
  mxf = prhs[4];

  if (!mxIsClass(prhs[5], "double"))
    mexErrMsgIdAndTxt("MCI:BadArgument", "calls must be nonnegative integer");
  calls = (int)mxGetScalar(prhs[5]);

  gslf.f = &integrand;
  gslf.dim = dim;
  gslf.params = &mxf; /* referencing because mxf is const */

  rng = gsl_rng_alloc(gsl_rng_default);

  if (strcmp(algorithm, "vegas") == 0) {
    gsl_monte_vegas_params params;
    gsl_monte_vegas_state *state;
    int i;
    double chisq;
            
    state = gsl_monte_vegas_alloc(dim);
    
    gsl_monte_vegas_params_get(state, &params);
    for (i = minnrhs; i + 1 < nrhs; i += 2) {
      char *k;
      double v;
      
      k = getNewString(prhs[i]);
      v = mxGetScalar(prhs[i+1]);
      
           if (strcmp(k, "alpha"     ) == 0) params.alpha      =         v;
      else if (strcmp(k, "iterations") == 0) params.iterations = (size_t)v;
      else if (strcmp(k, "stage"     ) == 0) params.stage      =    (int)v;
      else if (strcmp(k, "mode"      ) == 0) params.mode       =    (int)v;
      else if (strcmp(k, "verbose"   ) == 0) params.verbose    =    (int)v;
      else mexErrMsgIdAndTxt("MCI:BadArgument", "unknown parameter: %s", prhs[i]);
      
      free(k);
    }
    gsl_monte_vegas_params_set(state, &params);

    if (gsl_monte_vegas_integrate(&gslf, A, B, dim, calls, rng, state, &result, &abserr) == 0) {
      chisq = gsl_monte_vegas_chisq(state);
      if (abs(1 - chisq) > 1e3)
        mexWarnMsgIdAndTxt("MCI:ChisqInconsistent", "Chi-squared statistic is %f, which may be too far from 1.  Results may be inaccurate.", chisq);
    }

    gsl_monte_vegas_free(state);
  } else if (strcmp(algorithm, "miser") == 0) {
    gsl_monte_miser_params params;
    gsl_monte_miser_state *state;
    int i;

    state = gsl_monte_miser_alloc(dim);

    gsl_monte_miser_params_get(state, &params);
    for (i = minnrhs; i + 1 < nrhs; i += 2) {
      char *k;
      double v;
      
      k = getNewString(prhs[i]);
      v = mxGetScalar(prhs[i+1]);
      
           if (strcmp(k, "estimate_frac")           == 0) params.estimate_frac           =         v;
      else if (strcmp(k, "min_calls")               == 0) params.min_calls               = (size_t)v;
      else if (strcmp(k, "min_calls_per_bisection") == 0) params.min_calls_per_bisection = (size_t)v;
      else if (strcmp(k, "alpha")                   == 0) params.alpha                   =         v;
      else if (strcmp(k, "dither")                  == 0) params.dither                  =         v;
      else mexErrMsgIdAndTxt("MCI:BadArgument", "unknown parameter: %s", prhs[i]);
      
      free(k);
    }
    gsl_monte_miser_params_set(state, &params);

    gsl_monte_miser_integrate(&gslf, A, B, dim, calls, rng, state, &result, &abserr);

    gsl_monte_miser_free(state);
  } else if (strcmp(algorithm, "plain") == 0) {
    gsl_monte_plain_state *state = gsl_monte_plain_alloc(dim);
    gsl_monte_plain_integrate(&gslf, A, B, dim, calls, rng, state, &result, &abserr);
    gsl_monte_plain_free(state);
  } else {
    mexErrMsgIdAndTxt("MCI:BadArgument", "unknown algorithm: %s", algorithm);
  }
  
  free(algorithm);
  free(A);
  free(B);
  gsl_rng_free(rng);

  if (exception != NULL) {
    mxArray *ex = exception;
    exception = NULL;
    mexCallMATLAB(0, NULL, 1, &ex, "throw");
  }

  plhs[0] = mxCreateDoubleScalar(result);
  if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(abserr);
}
