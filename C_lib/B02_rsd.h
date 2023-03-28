#include <stdio.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

struct input{
 double *tr;
 double *cosm_par;
 double *pk_in;
 double sig_fog;
 gsl_interp_accel *acc_me;
 gsl_spline *spline_me;
 double *af;
 int mp;
};
