#include <stdio.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

struct input{
 double *tr;
 double *cosm_par;
 double *pk_in;
 double sig_fog;
 double sig_lin;
 gsl_interp_accel *acc;
 gsl_spline *spline;
 gsl_interp_accel *acc_n;
 gsl_spline *spline_n;
 gsl_interp_accel *acc_me;
 gsl_spline *spline_me;
 double knl;
 double D_z;
 int mod;
 double *af;
 double *ag;
 int mat;
 int mp;
};
struct input3{
 gsl_interp_accel *acc;
 gsl_spline *spline;
};
struct input2{
 double *tr;
 double *cosm_par;
 double sig_fog;
 double sig_lin;
 gsl_interp_accel *acc;
 gsl_spline *spline;
 double pka_l;
 double pkb_l;
 double pkc_l;
};

