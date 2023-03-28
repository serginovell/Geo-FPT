#include <stdio.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>


struct str_bkcbrt{
 double * kh_ar;
 double * par_ar;
 double * fit_ar;
 int idx_mup;
 gsl_interp_accel *acc;
 gsl_spline *spline;
};
 
