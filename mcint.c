#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_rng.h>

#include "mex.h"
#include "matrix.h"

char *getNewString(const mxArray *ms) {
  size_t len = mxGetNumberOfElements(ms);
  char *s = malloc(len * sizeof(char) + 1);
  mxGetString(ms, s, len + 1);
  return s;
}

double integrand(double x[], size_t dim, void *p) {
  mxArray *lhs[1], *rhs[2];
  rhs[0] = *((mxArray **)p);
  rhs[1] = mxCreateDoubleMatrix(1, dim, mxREAL);

  assert(mxIsClass(rhs[0], "function_handle"));

  /* We could do this without copying, but I'm not sure how safe it is to
   * change a Matlab object's data pointer to point to a region of memory
   * that's owned by the GSL routines.  It might be freed or the data might
   * be changed by GSL without Matlab knowing it. */
  double *x2 = mxGetPr(rhs[1]);
  for (int i = 0; i < dim; i++)
    x2[i] = x[i];

  mexCallMATLAB(1, lhs, 2, rhs, "feval");
  double y = mxGetScalar(lhs[0]);

  mxDestroyArray(rhs[1]);
  mxDestroyArray(lhs[0]);
  return y;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  const int minnrhs = 6;

  if (nrhs < minnrhs)
    mexErrMsgIdAndTxt("MCI:BadArgument", "not enough arguments given");

  if (!mxIsChar(prhs[0]))
    mexErrMsgIdAndTxt("MCI:BadArgument", "algorithm must be a string; use one of {plain,miser,vegas}");
  char *algorithm = getNewString(prhs[0]);

  if (!mxIsDouble(prhs[1]) || mxGetM(prhs[1]) != 1 || mxGetN(prhs[1]) != 1)
    mexErrMsgIdAndTxt("MCI:BadArgument", "dim must be a nonnegative integer");
  int dim = (int)mxGetScalar(prhs[1]);

  if (!mxIsDouble(prhs[2]) || mxGetM(prhs[2]) != 1 || mxGetN(prhs[2]) != dim)
    mexErrMsgIdAndTxt("MCI:BadArgument", "A must be a 1-by-dim vector specifying the lower bounds of each component of x");
  if (!mxIsDouble(prhs[3]) || mxGetM(prhs[3]) != 1 || mxGetN(prhs[3]) != dim)
    mexErrMsgIdAndTxt("MCI:BadArgument", "B must be a 1-by-dim vector specifying the upper bounds of each component of x");
  double *A = mxGetPr(prhs[2]), *B = mxGetPr(prhs[3]);

  if (!mxIsClass(prhs[4], "function_handle"))
    mexErrMsgIdAndTxt("MCI:BadArgument", "f must be a function handle");
  const mxArray *mxf = prhs[4];

  if (!mxIsClass(prhs[5], "double"))
    mexErrMsgIdAndTxt("MCI:BadArgument", "calls must be nonnegative integer");
  int calls = (int)mxGetScalar(prhs[5]);

  gsl_monte_function gslf;
  gslf.f = &integrand;
  gslf.dim = dim;
  gslf.params = &mxf; /* referencing because mxf is const */

  gsl_rng *rng = gsl_rng_alloc(gsl_rng_default);

  double result = 0, abserr = 0;

  if (strcmp(algorithm, "vegas") == 0) {
    gsl_monte_vegas_state *state = gsl_monte_vegas_alloc(dim);
    gsl_monte_vegas_params params;

    gsl_monte_vegas_params_get(state, &params);
    for (int i = minnrhs; i + 1 < nrhs; i += 2) {
      char *k = getNewString(prhs[i]);
      double v = mxGetScalar(prhs[i+1]);
           if (strcmp(k, "alpha"     ) == 0) params.alpha      =         v;
      else if (strcmp(k, "iterations") == 0) params.iterations = (size_t)v;
      else if (strcmp(k, "stage"     ) == 0) params.stage      =    (int)v;
      else if (strcmp(k, "mode"      ) == 0) params.mode       =    (int)v;
      else if (strcmp(k, "verbose"   ) == 0) params.verbose    =    (int)v;
      else mexErrMsgIdAndTxt("MCI:BadArgument", "unknown parameter: %s", prhs[i]);
    }
    gsl_monte_vegas_params_set(state, &params);

    gsl_monte_vegas_integrate(&gslf, A, B, dim, calls, rng, state, &result, &abserr);

    double chisq = gsl_monte_vegas_chisq(state);
    if (abs(1 - chisq) > 1e3)
      mexWarnMsgIdAndTxt("MCI:ChisqInconsistent", "Chi-squared statistic is %f, which may be too far from 1.  Results may be inaccurate.", chisq);

    gsl_monte_vegas_free(state);
  } else if (strcmp(algorithm, "miser") == 0) {
    gsl_monte_miser_state *state = gsl_monte_miser_alloc(dim);
    gsl_monte_miser_params params;

    gsl_monte_miser_params_get(state, &params);
    for (int i = minnrhs; i + 1 < nrhs; i += 2) {
      char *k = getNewString(prhs[i]);
      double v = mxGetScalar(prhs[i+1]);
           if (strcmp(k, "estimate_frac")           == 0) params.estimate_frac           =         v;
      else if (strcmp(k, "min_calls")               == 0) params.min_calls               = (size_t)v;
      else if (strcmp(k, "min_calls_per_bisection") == 0) params.min_calls_per_bisection = (size_t)v;
      else if (strcmp(k, "alpha")                   == 0) params.alpha                   =         v;
      else if (strcmp(k, "dither")                  == 0) params.dither                  =         v;
      else mexErrMsgIdAndTxt("MCI:BadArgument", "unknown parameter: %s", prhs[i]);
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

  if (nlhs > 0) plhs[0] = mxCreateDoubleScalar(result);
  if (nlhs > 1) plhs[1] = mxCreateDoubleScalar(abserr);
}