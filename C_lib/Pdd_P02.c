#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <complex.h>//complex.h always BEFORE fftw3.h
#include <fftw3.h>
#include <omp.h>

#include "functions.h"
#include "cubature.h"
#include "Pdd_P02.h"

double N2(double pl,double p13,double p15)
{
        double f;
    if (p15>=0)
    {
                f=cosh(sqrt(p15/pl))+0.5*p13/pl*sqrt(pl/p15)*sinh(sqrt(p15/pl));
    }
        else
        {
        f=cos(sqrt(-p15/pl))+0.5*p13/pl*sqrt(-pl/p15)*sin(sqrt(-p15/pl));
        }

        f=f*f;
        return f;
}


void Preal(double *kin,double *Pin,int Nin, void *fdata)
{

        f_params params_function = *(f_params *) fdata;

    int i,N;
    double k1p,aiso,Pdd,alpha_perpendicular,alpha_parallel,sigma8_scaling,mBGV,m2BGV;
    double **Theory;
    char *ptmodel;
    int N1;
    double w1,w0,w2;
    int interpolation_order,shiftN;
    char *spacing;
    double shape_factor,s8;
    double kpivot=0.03;//this scale corresponds to 8Mpc/h
    double ascale=0.6;//scale at which the asymptotes are reached


    interpolation_order=1;//1 or 2
    if(interpolation_order==1){shiftN=1;}
    if(interpolation_order==2){shiftN=2;}

    spacing=params_function.spacing;
    N=params_function.N;
    alpha_perpendicular=params_function.a_perpendicular;
    alpha_parallel=params_function.a_parallel;
    mBGV=params_function.mBGV;
    m2BGV=params_function.m2BGV;
    sigma8_scaling=params_function.sigma8;
    Theory=params_function.theory;
    ptmodel=params_function.type_ptmodel;

aiso=pow(alpha_perpendicular*alpha_perpendicular*alpha_parallel,1./3.);
for(i=0;i<Nin;i++)
{

k1p=kin[i];
//shape_factor=1+mBGV*(log10(k1p)-log10(kpivot));
shape_factor=exp(  mBGV*log(k1p/kpivot)+m2BGV/ascale*tanh(ascale*log(k1p/kpivot)) );
s8=params_function.sigma8;
sigma8_scaling=shape_factor*s8;//trial option. Put shape factor at the same level as s8**2

N1=determine_N_doublearray(Theory,k1p,N,spacing);
if(N1>=N-shiftN || N1<0 || k1p<=0 || kin[i]<=0){Pdd=0;}
else{

if(interpolation_order==1){
w1=determine_w1_doublearray(Theory,k1p,N1,spacing);
}
if(interpolation_order==2){
w0=determine_w0_2ndorder_doublearray(Theory,k1p,N1,spacing);
w1=determine_w1_2ndorder_doublearray(Theory,k1p,N1,spacing);
w2=determine_w2_2ndorder_doublearray(Theory,k1p,N1,spacing);
}

if( strcmp(ptmodel, "linear") == 0){Pdd = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2);}

if( strcmp(ptmodel, "1L-SPT") == 0){Pdd = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2));}

if( strcmp(ptmodel, "2L-SPT") == 0){Pdd = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2))+sigma8_scaling*sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,5,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,12,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,10,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)/(4.*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)));}

if( strcmp(ptmodel, "1L-RPT") == 0){Pdd = (sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2))*exp(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)/P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2));}

if( strcmp(ptmodel, "2L-RPT") == 0){Pdd =(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,5,N,spacing,interpolation_order,N1,w0,w1,w2))*N2(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,10,N,spacing,interpolation_order,N1,w0,w1,w2));}

}

Pin[i]=Pdd;

}

}



int integralP(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{

    struct f_params2 * fp = (struct f_params2 *)fdata;
    
    
    int N, nomask,noise_option;
    double b1,b2,bs2,b3nl;//,Anoise,Pnoise;
    double alpha_perpendicular,alpha_parallel,F,sigma8_scaling,sigma_ps,mBGV,m2BGV,s8,kpivot,shape_factor,ascale;
    double **Theory;
    double a11,a12,a22,a23,a33;
    double b1_11,b1_12,b1_21,b1_22;
    double b2_11,b2_12,b2_21,b2_22;
    double b3_12,b3_21,b3_22,b4_22;
    double Legendre, Pdd, Ptt, Pdt, Pb2_d, Pb2_t, Pbs2_d, Pbs2_t, sigma3, Pb22, Pb2s2, Pbs22;
    double fog;
    char *fog_model_ps, *ptmodel, *rsdmodel_ps;
    int mode;
    double kinput, mu, Fap, mup;
    double k1p;
    int N1;
    double w1,w0,w2;
    int interpolation_order,shiftN;
    char *spacing;
    kpivot=0.03;//8Mpc/h
    ascale=0.6;
    //remove warning
    (void)ndim;
    (void)fdim;
    interpolation_order=1;//1 or 2
    if(interpolation_order==1){shiftN=1;}
    if(interpolation_order==2){shiftN=2;}

    spacing= fp-> spacing;
    N= fp-> N;
    b1= fp-> b1;
    b2= fp-> b2;
    bs2= fp-> bs2;
    b3nl= fp-> b3nl;
    //Anoise= fp-> A;
    //Pnoise= fp-> Pnoise;
    alpha_perpendicular= fp-> a_perpendicular;
    alpha_parallel= fp-> a_parallel;
    F= fp-> f;
    mBGV= fp-> mBGV;
    m2BGV= fp-> m2BGV;
//    sigma8_scaling=params_function.sigma8;
    s8= fp-> sigma8;
    sigma_ps= fp-> sigmaP;
    Theory= fp-> theory;

    fog_model_ps= fp-> type_fog;
    ptmodel= fp-> type_ptmodel;
    rsdmodel_ps= fp-> type_rsdmodel;
    noise_option=0;//params_function.noise_option;
    mode= fp-> mode;
    kinput= fp-> kin;
    mu=x[0];
    Fap=alpha_parallel*1./alpha_perpendicular*1.;
    mup=mu/Fap*pow(1.+mu*mu*(1./(Fap*Fap)-1),-0.5);
    k1p=kinput/alpha_perpendicular*pow( 1+mu*mu*(1./(Fap*Fap)-1) ,0.5);

    N1=determine_N_doublearray(Theory,k1p,N,spacing);

 //   shape_factor=1+mBGV*(log10(k1p)-log10(kpivot));
      shape_factor=exp(  mBGV*log(k1p/kpivot)+m2BGV/ascale*tanh(ascale*log(k1p/kpivot)) );
  sigma8_scaling=shape_factor*s8;//trial option. Put shape factor at the same level as s8**2

if(N1>=N-shiftN || N1<0 || k1p<=0 || kinput<=0)
{
fval[0]=0;//integral returns 0 out of range of Theory
}
else{
    //w1=determine_w1(Theory,k1p,N1);

if(interpolation_order==1){
    w1=determine_w1_doublearray(Theory,k1p,N1,spacing);
}
if(interpolation_order==2){

w0=determine_w0_2ndorder_doublearray(Theory,k1p,N1,spacing);
w1=determine_w1_2ndorder_doublearray(Theory,k1p,N1,spacing);
w2=determine_w2_2ndorder_doublearray(Theory,k1p,N1,spacing);


}

//if( strcmp(rsdmodel_ps, "TNS10") ==0){

    a11=P_interpol_fast_doublearray(k1p,Theory,24,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*b1*F;
    a12=P_interpol_fast_doublearray(k1p,Theory,25,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F;
    a22=P_interpol_fast_doublearray(k1p,Theory,26,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F;
    a23=P_interpol_fast_doublearray(k1p,Theory,27,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F;
    a33=P_interpol_fast_doublearray(k1p,Theory,28,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F;

    b1_11=P_interpol_fast_doublearray(k1p,Theory,29,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*b1*F*F;
    b1_12=P_interpol_fast_doublearray(k1p,Theory,30,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b1_21=P_interpol_fast_doublearray(k1p,Theory,31,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b1_22=P_interpol_fast_doublearray(k1p,Theory,32,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F*F;

    b2_11=P_interpol_fast_doublearray(k1p,Theory,33,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*b1*F*F;
    b2_12=P_interpol_fast_doublearray(k1p,Theory,34,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b2_21=P_interpol_fast_doublearray(k1p,Theory,35,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b2_22=P_interpol_fast_doublearray(k1p,Theory,36,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F*F;

    b3_12=P_interpol_fast_doublearray(k1p,Theory,37,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b3_21=P_interpol_fast_doublearray(k1p,Theory,38,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*b1*F*F*F;
    b3_22=P_interpol_fast_doublearray(k1p,Theory,39,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F*F;
    b4_22=P_interpol_fast_doublearray(k1p,Theory,40,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling*F*F*F*F;

//}


if( strcmp(ptmodel, "linear") == 0)
{

Pdd = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2);
Pdt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2);
Ptt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2);

}

if( strcmp(ptmodel, "1L-SPT") == 0)
{

Pdd = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2));
Pdt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,3,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2));
Ptt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,4,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2));

}


if( strcmp(ptmodel, "2L-SPT") == 0)
{

Pdd = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2))+sigma8_scaling*sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,5,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,12,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,10,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)/(4.*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)));

Pdt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,3,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2))+sigma8_scaling*sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,6,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,14,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,15,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,10,N,spacing,interpolation_order,N1,w0,w1,w2)+0.5*P_interpol_fast_doublearray(k1p,Theory,11,N,spacing,interpolation_order,N1,w0,w1,w2)+0.25*(P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2))*(P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2))/(4.*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)));

Ptt = sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,4,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2))+sigma8_scaling*sigma8_scaling*sigma8_scaling*(P_interpol_fast_doublearray(k1p,Theory,7,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,13,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,11,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2)*P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2)/(4.*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)));


}

if( strcmp(ptmodel, "1L-RPT") == 0)
{

Pdd = (sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2))*exp(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)/P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2));
Pdt = (sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,3,N,spacing,interpolation_order,N1,w0,w1,w2))*exp(sigma8_scaling*0.5*(P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2))/P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2));
Ptt = (sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,4,N,spacing,interpolation_order,N1,w0,w1,w2))*exp(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2)/P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2));

}


if( strcmp(ptmodel, "2L-RPT") == 0)
{

Pdd =(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,2,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,5,N,spacing,interpolation_order,N1,w0,w1,w2))*N2(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,10,N,spacing,interpolation_order,N1,w0,w1,w2));
Pdt = (sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,3,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,6,N,spacing,interpolation_order,N1,w0,w1,w2))*N2(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*0.5*(P_interpol_fast_doublearray(k1p,Theory,8,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2)),sigma8_scaling*sigma8_scaling*sigma8_scaling*0.5*(P_interpol_fast_doublearray(k1p,Theory,10,N,spacing,interpolation_order,N1,w0,w1,w2)+P_interpol_fast_doublearray(k1p,Theory,11,N,spacing,interpolation_order,N1,w0,w1,w2)));
Ptt = (sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,4,N,spacing,interpolation_order,N1,w0,w1,w2)+sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,7,N,spacing,interpolation_order,N1,w0,w1,w2))*N2(sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,1,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,9,N,spacing,interpolation_order,N1,w0,w1,w2),sigma8_scaling*sigma8_scaling*sigma8_scaling*P_interpol_fast_doublearray(k1p,Theory,11,N,spacing,interpolation_order,N1,w0,w1,w2));


}



     Pb2_d=P_interpol_fast_doublearray(k1p,Theory,16,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pb2_t=P_interpol_fast_doublearray(k1p,Theory,17,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pbs2_d=P_interpol_fast_doublearray(k1p,Theory,18,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pbs2_t=P_interpol_fast_doublearray(k1p,Theory,19,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     sigma3=P_interpol_fast_doublearray(k1p,Theory,23,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pb22=P_interpol_fast_doublearray(k1p,Theory,22,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pb2s2=P_interpol_fast_doublearray(k1p,Theory,20,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;
     Pbs22=P_interpol_fast_doublearray(k1p,Theory,21,N,spacing,interpolation_order,N1,w0,w1,w2)*sigma8_scaling*sigma8_scaling;


if( strcmp(fog_model_ps, "Lorentzian") ==0){fog=pow( 1+sigma_ps*sigma_ps/2.*k1p*k1p*mup*mup ,-2);}
if( strcmp(fog_model_ps, "Exponential") ==0){fog=exp( -sigma_ps*sigma_ps/2.*k1p*k1p*mup*mup);}

if(fabs(sigma_ps)<0.001){fog=1.;}

if(mode==0 && noise_option == 0){Legendre=1.0;nomask=1;}
if(mode==0 && noise_option == 1 ){Legendre=1.0;nomask=0;}
if(mode==2){Legendre=3./2.*mu*mu-0.5;nomask=0;}
if(mode==4){Legendre=35./8.*mu*mu*mu*mu-30./8.*mu*mu+3./8.;nomask=0;}


if( strcmp(rsdmodel_ps, "TNS10") ==0){

              fval[0]=(2*mode+1)/2.*((b1*b1*Pdd+2.*b2*b1*Pb2_d+2.*b1*bs2*Pbs2_d+2.*sigma3*b1*b3nl+b2*b2*Pb22+2.*b2*bs2*Pb2s2+bs2*bs2*Pbs22)+2.*F*mup*mup*( b1*Pdt+b2*Pb2_t+bs2*Pbs2_t+sigma3*b3nl )+pow(mup,4)*F*F*Ptt+a11*mup*mup+a12*mup*mup+a22*mup*mup*mup*mup+a23*mup*mup*mup*mup+a33*mup*mup*mup*mup*mup*mup+b1_11*mup*mup+b1_12*mup*mup+b1_21*mup*mup+b1_22*mup*mup+b2_11*mup*mup*mup*mup+b2_12*mup*mup*mup*mup+b2_21*mup*mup*mup*mup+b2_22*mup*mup*mup*mup+b3_12*mup*mup*mup*mup*mup*mup+b3_21*mup*mup*mup*mup*mup*mup+b3_22*mup*mup*mup*mup*mup*mup+b4_22*mup*mup*mup*mup*mup*mup*mup*mup)*fog*Legendre/(alpha_parallel*alpha_perpendicular*alpha_perpendicular); //+nomask*Anoise*Pnoise*(2*mode+1)/(2*alpha_parallel*alpha_perpendicular*alpha_perpendicular)*Legendre;

}


if( strcmp(rsdmodel_ps, "Scoccimarro04") ==0){
              fval[0]=(2*mode+1)/2.*((b1*b1*Pdd+2.*b2*b1*Pb2_d+2.*b1*bs2*Pbs2_d+2.*sigma3*b1*b3nl+b2*b2*Pb22+2.*b2*bs2*Pb2s2+bs2*bs2*Pbs22)+2.*F*mup*mup*( b1*Pdt+b2*Pb2_t+bs2*Pbs2_t+sigma3*b3nl )+pow(mup,4)*F*F*Ptt)*fog*Legendre/(alpha_parallel*alpha_perpendicular*alpha_perpendicular);//+nomask*Anoise*Pnoise*(2*mode+1)/(2*alpha_parallel*alpha_perpendicular*alpha_perpendicular)*Legendre;

}

if( strcmp(rsdmodel_ps, "Kaiser87") ==0){

 fval[0]=(2*mode+1)/2.*((b1*b1*Pdd+2.*b2*b1*Pb2_d+2.*b1*bs2*Pbs2_d+2.*sigma3*b1*b3nl+b2*b2*Pb22+2.*b2*bs2*Pb2s2+bs2*bs2*Pbs22)*(1+F/b1*mup*mup)*(1+F/b1*mup*mup))*fog*Legendre/(alpha_parallel*alpha_perpendicular*alpha_perpendicular);//+nomask*Anoise*Pnoise*(2*mode+1)/(2*alpha_parallel*alpha_perpendicular*alpha_perpendicular)*Legendre;

}

}
return 0;
}

extern void ext_Prsd(double **Theory, double *kin, int kin_dim, double *Pout, double alpa, double alpe, double sigma8_scaling, int Num,double b1,double b2,double bs2,double b3nl,double f,int mup,double sigmaP)
{
    double val,err;
    int i;
    double xmin[1]={-1};
    double xmax[1]={+1};
    for (i=0; i<kin_dim; i++)
    {
        struct f_params2 f_par = {"2L-RPT",Num,0.,0.,sigma8_scaling,sigma8_scaling,alpa,alpe,Theory,"linear",
            b1,b2,bs2,b3nl,f,"Lorentzian","TNS10",mup,kin[i],sigmaP};
        
        hcubature(1,integralP,&f_par,1,xmin,xmax,0, 0, 1e-5, 1e-5, &val, &err);
        Pout[i] = val;
    }    

}

  
 
extern void ext_Preal(double **Theory, double *kin, double *Pout, double alpa, double alpe, double sigma8_scaling, int Num)
{
                          
  
    f_params *f_par;
    f_par= (f_params *) malloc(sizeof(f_params));
    (*f_par).type_ptmodel="2L-RPT";
    (*f_par).N=Num;
    (*f_par).mBGV=0.;
    (*f_par).m2BGV=0.;
    (*f_par).sigma8=sigma8_scaling;
//    (*f_par).sigma8_value=sigma8_scaling;
    (*f_par).a_parallel=alpa;
    (*f_par).a_perpendicular=alpe;
    (*f_par).theory=Theory;
    (*f_par).spacing="linear";
    
    
    Preal(kin,Pout,Num,f_par);
}
