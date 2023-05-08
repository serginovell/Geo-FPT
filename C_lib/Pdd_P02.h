#include <stdio.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

void Preal(double *kin,double *Pin,int N, void *fdata);
extern void ext_Preal(double **Theory, double *kin, double *Pout, double alpa, double alpe, double sigma8_scaling, int Num);

typedef struct f_params {


char *type_ptmodel;
int N;
double mBGV;
double m2BGV;
double sigma8;
//double sigma8_value;
double a_parallel;
double a_perpendicular;
double **theory;
char *spacing;

}f_params;

typedef struct f_params2{
char *type_ptmodel;
int N;
double mBGV;
double m2BGV;
double sigma8;
double sigma8_value;
double a_parallel;
double a_perpendicular;
double **theory;
char *spacing;
double b1;
double b2;
double bs2;
double b3nl;
//double A;
//double Pnoise;
double f;
char *type_fog;
char *type_rsdmodel;
int mode;
double kin;
double sigmaP;
}f_params2;

