#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include "cubature.h"

typedef struct f_params3{
    double alpa;
    double alpe;
    double s8_sc;
    double b1;
    double b2;
    double F;
    double sP;
    double kknl;
    int multip;
    gsl_interp_accel *acc1;
    gsl_spline *spline1;
    gsl_interp_accel *acc2;
    gsl_spline *spline2;
    gsl_interp_accel *acc3;
    gsl_spline *spline3;
    gsl_interp_accel *acc4;
    gsl_spline *spline4;
    gsl_interp_accel *acc5;
    gsl_spline *spline5;
    gsl_interp_accel *acc6;
    gsl_spline *spline6;
    gsl_interp_accel *acc7;
    gsl_spline *spline7;
    gsl_interp_accel *acc8;
    gsl_spline *spline8;
    gsl_interp_accel *acc9;
    gsl_spline *spline9;
    gsl_interp_accel *acc10;
    gsl_spline *spline10;
    gsl_interp_accel *acc11;
    gsl_spline *spline11;
    gsl_interp_accel *acc12;
    gsl_spline *spline12;
    gsl_interp_accel *acc13;
    gsl_spline *spline13;
    gsl_interp_accel *acc14;
    gsl_spline *spline14;
    gsl_interp_accel *acc15;
    gsl_spline *spline15;
    gsl_interp_accel *acc16;
    gsl_spline *spline16;
    gsl_interp_accel *acc17;
    gsl_spline *spline17;
    gsl_interp_accel *acc18;
    gsl_spline *spline18;
    gsl_interp_accel *acc19;
    gsl_spline *spline19;
    gsl_interp_accel *acc20;
    gsl_spline *spline20;
    gsl_interp_accel *acc21;
    gsl_spline *spline21;
    gsl_interp_accel *acc22;
    gsl_spline *spline22;
    gsl_interp_accel *acc23;
    gsl_spline *spline23;
    gsl_interp_accel *acc24;
    gsl_spline *spline24;
    gsl_interp_accel *acc25;
    gsl_spline *spline25;
    gsl_interp_accel *acc26;
    gsl_spline *spline26;
    gsl_interp_accel *acc27;
    gsl_spline *spline27;
    gsl_interp_accel *acc28;
    gsl_spline *spline28;
    gsl_interp_accel *acc29;
    gsl_spline *spline29;
    gsl_interp_accel *acc30;
    gsl_spline *spline30;
    gsl_interp_accel *acc31;
    gsl_spline *spline31;
    gsl_interp_accel *acc32;
    gsl_spline *spline32;
    gsl_interp_accel *acc33;
    gsl_spline *spline33;
    gsl_interp_accel *acc34;
    gsl_spline *spline34;
    gsl_interp_accel *acc35;
    gsl_spline *spline35;
    gsl_interp_accel *acc36;
    gsl_spline *spline36;
    gsl_interp_accel *acc37;
    gsl_spline *spline37;
    gsl_interp_accel *acc38;
    gsl_spline *spline38;
    gsl_interp_accel *acc39;
    gsl_spline *spline39;
    gsl_interp_accel *acc40;
    gsl_spline *spline40;
    
}f_params3;
 
