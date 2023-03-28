#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include "cubature.h"
#include "Preal2.h"

int P024_int(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval)
{
    struct f_params3 * fp = (struct f_params3 *)fdata;

    double alpa,alpe,s8_sc,b1,b2,F,sP,kknl,mu2,Leg;
    int multip;
    (void)ndim;
    (void)fdim;
    alpa = fp->alpa;
    alpe = fp->alpe;
    s8_sc = fp->s8_sc;
    b1 = fp->b1;
    b2 = fp->b2;
    F = fp->F;
    sP = fp->sP;
    mu2 = x[0];
    kknl = fp->kknl;
    multip = fp-> multip;
    gsl_interp_accel *acc1 = fp->acc1;
    gsl_spline *spline1    = fp->spline1;
    gsl_interp_accel *acc2 = fp->acc2;
    gsl_spline *spline2    = fp->spline2;
    gsl_interp_accel *acc3 = fp->acc3;
    gsl_spline *spline3    = fp->spline3;
    gsl_interp_accel *acc4 = fp->acc4;
    gsl_spline *spline4    = fp->spline4;
    gsl_interp_accel *acc5 = fp->acc5;
    gsl_spline *spline5    = fp->spline5;
    gsl_interp_accel *acc6 = fp->acc6;
    gsl_spline *spline6    = fp->spline6;
    gsl_interp_accel *acc7 = fp->acc7;
    gsl_spline *spline7    = fp->spline7;
    gsl_interp_accel *acc8 = fp->acc8;
    gsl_spline *spline8    = fp->spline8;
    gsl_interp_accel *acc9 = fp->acc9;
    gsl_spline *spline9    = fp->spline9;
    gsl_interp_accel *acc10 = fp->acc10;
    gsl_spline *spline10    = fp->spline10;
    gsl_interp_accel *acc11 = fp->acc11;
    gsl_spline *spline11    = fp->spline11;
    gsl_interp_accel *acc12 = fp->acc12;
    gsl_spline *spline12    = fp->spline12;
    gsl_interp_accel *acc13 = fp->acc13;
    gsl_spline *spline13    = fp->spline13;
    gsl_interp_accel *acc14 = fp->acc14;
    gsl_spline *spline14    = fp->spline14;
    gsl_interp_accel *acc15 = fp->acc15;
    gsl_spline *spline15    = fp->spline15;
    gsl_interp_accel *acc16 = fp->acc16;
    gsl_spline *spline16    = fp->spline16;
    gsl_interp_accel *acc17 = fp->acc17;
    gsl_spline *spline17    = fp->spline17;
    gsl_interp_accel *acc18 = fp->acc18;
    gsl_spline *spline18    = fp->spline18;
    gsl_interp_accel *acc19 = fp->acc19;
    gsl_spline *spline19    = fp->spline19;
    gsl_interp_accel *acc20 = fp->acc20;
    gsl_spline *spline20    = fp->spline20;
    gsl_interp_accel *acc21 = fp->acc21;
    gsl_spline *spline21    = fp->spline21;
    gsl_interp_accel *acc22 = fp->acc22;
    gsl_spline *spline22    = fp->spline22;
    gsl_interp_accel *acc23 = fp->acc23;
    gsl_spline *spline23    = fp->spline23;
    gsl_interp_accel *acc24 = fp->acc24;
    gsl_spline *spline24    = fp->spline24;
    gsl_interp_accel *acc25 = fp->acc25;
    gsl_spline *spline25    = fp->spline25;
    gsl_interp_accel *acc26 = fp->acc26;
    gsl_spline *spline26    = fp->spline26;
    gsl_interp_accel *acc27 = fp->acc27;
    gsl_spline *spline27    = fp->spline27;
    gsl_interp_accel *acc28 = fp->acc28;
    gsl_spline *spline28    = fp->spline28;
    gsl_interp_accel *acc29 = fp->acc29;
    gsl_spline *spline29    = fp->spline29;
    gsl_interp_accel *acc30 = fp->acc30;
    gsl_spline *spline30    = fp->spline30;
    gsl_interp_accel *acc31 = fp->acc31;
    gsl_spline *spline31    = fp->spline31;
    gsl_interp_accel *acc32 = fp->acc32;
    gsl_spline *spline32    = fp->spline32;
    gsl_interp_accel *acc33 = fp->acc33;
    gsl_spline *spline33    = fp->spline33;
    gsl_interp_accel *acc34 = fp->acc34;
    gsl_spline *spline34    = fp->spline34;
    gsl_interp_accel *acc35 = fp->acc35;
    gsl_spline *spline35    = fp->spline35;
    gsl_interp_accel *acc36 = fp->acc36;
    gsl_spline *spline36    = fp->spline36;
    gsl_interp_accel *acc37 = fp->acc37;
    gsl_spline *spline37    = fp->spline37;
    gsl_interp_accel *acc38 = fp->acc38;
    gsl_spline *spline38    = fp->spline38;
    gsl_interp_accel *acc39 = fp->acc39;
    gsl_spline *spline39    = fp->spline39;
    gsl_interp_accel *acc40 = fp->acc40;
    gsl_spline *spline40    = fp->spline40;
//    gsl_interp_accel *acc41 = fp->acc41;
    //gsl_spline *spline41    = fp->spline41;
    
    
    double bs2 = -4./7.*(b1-1.);
    double b3nl= 32./315.*(b1-1.);
    double Fap=alpa/alpe;
    double mup=mu2/(Fap*sqrt(1.+mu2*mu2*(1./(Fap*Fap)-1.)));
    double k1p=kknl/alpe*sqrt(1.+mu2*mu2*(1./(Fap*Fap)-1));
    double logk = log10(k1p);
    double P1= gsl_spline_eval(spline1, logk, acc1);
    double P2= gsl_spline_eval(spline2, logk, acc2);
    double P3= gsl_spline_eval(spline3, logk, acc3);
    double P4= gsl_spline_eval(spline4, logk, acc4);
    double P5= gsl_spline_eval(spline5, logk, acc5);
    double P6= gsl_spline_eval(spline6, logk, acc6);
    double P7= gsl_spline_eval(spline7, logk, acc7);
    double P8= gsl_spline_eval(spline8, logk, acc8);
    double P9= gsl_spline_eval(spline9, logk, acc9);
    double P10= gsl_spline_eval(spline10, logk, acc10);
    double P11= gsl_spline_eval(spline11, logk, acc11);
    double P12= gsl_spline_eval(spline12, logk, acc12);
    double P13= gsl_spline_eval(spline13, logk, acc13);
    double P14= gsl_spline_eval(spline14, logk, acc14);
    double P15= gsl_spline_eval(spline15, logk, acc15);
    double P16= gsl_spline_eval(spline16, logk, acc16);
    double P17= gsl_spline_eval(spline17, logk, acc17);
    double P18= gsl_spline_eval(spline18, logk, acc18);
    double P19= gsl_spline_eval(spline19, logk, acc19);
    double P20= gsl_spline_eval(spline20, logk, acc20);
    double P21= gsl_spline_eval(spline21, logk, acc21);
    double P22= gsl_spline_eval(spline22, logk, acc22);
    double P23= gsl_spline_eval(spline23, logk, acc23);
    double P24= gsl_spline_eval(spline24, logk, acc24);
    double P25= gsl_spline_eval(spline25, logk, acc25);
    double P26= gsl_spline_eval(spline26, logk, acc26);
    double P27= gsl_spline_eval(spline27, logk, acc27);
    double P28= gsl_spline_eval(spline28, logk, acc28);
    double P29= gsl_spline_eval(spline29, logk, acc29);
    double P30= gsl_spline_eval(spline30, logk, acc30);
    double P31= gsl_spline_eval(spline31, logk, acc31);
    double P32= gsl_spline_eval(spline32, logk, acc32);
    double P33= gsl_spline_eval(spline33, logk, acc33);
    double P34= gsl_spline_eval(spline34, logk, acc34);
    double P35= gsl_spline_eval(spline35, logk, acc35);
    double P36= gsl_spline_eval(spline36, logk, acc36);
    double P37= gsl_spline_eval(spline37, logk, acc37);
    double P38= gsl_spline_eval(spline38, logk, acc38);
    double P39= gsl_spline_eval(spline39, logk, acc39);
    double P40= gsl_spline_eval(spline40, logk, acc40);
    //double P41= gsl_spline_eval(spline41, logk, acc41);
    

    double a11=P24*s8_sc*s8_sc*b1*b1*F,a12=P25*s8_sc*s8_sc*b1*F*F,a22=P26*s8_sc*s8_sc*b1*F*F,a23=P27*s8_sc*s8_sc*F*F*F,a33 = P28*F*F*F;
    double b1_11=P29*s8_sc*s8_sc*b1*b1*F*F,b1_12=P30*s8_sc*s8_sc*b1*F*F*F,b1_21=P31*s8_sc*s8_sc*b1*F*F*F,b1_22=P32*s8_sc*s8_sc*F*F*F*F;
    double b2_11=P33*s8_sc*s8_sc*b1*b1*F*F,b2_12=P34*s8_sc*s8_sc*b1*F*F*F,b2_21=P35*s8_sc*s8_sc*b1*F*F*F,b2_22=P36*s8_sc*s8_sc*F*F*F*F;
    double b3_12=P37*s8_sc*s8_sc*b1*F*F*F,b3_21=P38*s8_sc*s8_sc*b1*F*F*F,b3_22=P39*s8_sc*s8_sc*F*F*F*F,b4_22=P40*s8_sc*s8_sc*F*F*F*F;

    
    double Pddd   = s8_sc*P1+s8_sc*s8_sc*(P2+P8)+pow(s8_sc,3)*(P5+P12+P10+P8*P8/(4.*P1));
    double Pdt    = s8_sc*P1+s8_sc*s8_sc*(P3+0.5*(P8+P9))+pow(s8_sc,3)*(P6+0.5*(P14+P15+P10+P11)+0.25*pow(P8+P9,2)/(4.*P1));
    double Ptt    = s8_sc*P1+s8_sc*s8_sc*(P4+P9)+pow(s8_sc,3)*(P7+P13+P11+P9*P9/(4.*P1));    
    double Pb2_d  = pow(s8_sc,2)*P16;
    double Pb2_t  = pow(s8_sc,2)*P17;
    double Pbs2_d = pow(s8_sc,2)*P18;
    double Pbs2_t = pow(s8_sc,2)*P19;
    double sigma3 = pow(s8_sc,2)*P23;
    double Pb22   = pow(s8_sc,2)*P22;
    double Pb2s2  = pow(s8_sc,2)*P20;
    double Pbs22  = pow(s8_sc,2)*P21;
    
    double fog = 1./pow( 1+0.5*pow(sP*k1p*mup,2),2 );
    if(multip==0){ Leg = 1.;}
    if(multip==2){Leg = 0.5*(3.*mu2*mu2-1.);}
    if(multip==4){Leg = 1./8.*(35.*pow(mu2,4)-30.*mu2*mu2+3);}
        

    fval[0] = (2*multip+1)/2.*((b1*b1*Pddd+2.*b2*b1*Pb2_d+2.*b1*bs2*Pbs2_d+2.*sigma3*b1*b3nl+b2*b2*Pb22+2.*b2*bs2*Pb2s2+bs2*bs2*Pbs22)+2.*F*mup*mup*( b1*Pdt+b2*Pb2_t+bs2*Pbs2_t+sigma3*b3nl )+pow(mup,4)*F*F*Ptt+a11*mup*mup+a12*mup*mup+a22*mup*mup*mup*mup+a23*mup*mup*mup*mup+a33*mup*mup*mup*mup*mup*mup+b1_11*mup*mup+b1_12*mup*mup+b1_21*mup*mup+b1_22*mup*mup+b2_11*mup*mup*mup*mup+b2_12*mup*mup*mup*mup+b2_21*mup*mup*mup*mup+b2_22*mup*mup*mup*mup+b3_12*mup*mup*mup*mup*mup*mup+b3_21*mup*mup*mup*mup*mup*mup+b3_22*mup*mup*mup*mup*mup*mup+b4_22*mup*mup*mup*mup*mup*mup*mup*mup)*fog*Leg/(alpa*alpe*alpe);
    
    return 0; 
}


extern void ext_P024(double **Theory, double *kknl, int kknl_dim, int th_dim, double *Pout, double *cosm_par,int multip,double *k_t,double *P1,double *P2,double *P3,double *P4,double *P5,double *P6,double *P7,double *P8,double *P9,double *P10, double *P11,double *P12,double *P13,double *P14,double *P15,double *P16,double *P17,double *P18,double *P19,double *P20,double *P21,double *P22,double *P23,double *P24,double *P25,double *P26,double *P27,double *P28,double *P29,double *P30, double *P31,double *P32,double *P33,double *P34,double *P35,double *P36,double *P37,double *P38,double *P39,double *P40)
{
    
    double val,err;
    double alpa=cosm_par[5],alpe=cosm_par[6],s8_sc=pow(cosm_par[0]/Theory[0][41],2),b1=cosm_par[1],b2=cosm_par[3],F=cosm_par[2],sP=cosm_par[9];
    int i;
    double xmin[1]={-1};
    double xmax[1]={+1};

    gsl_interp_accel *acc1 = gsl_interp_accel_alloc ();
    gsl_spline *spline1 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline1, k_t, P1, th_dim);
    gsl_interp_accel *acc2 = gsl_interp_accel_alloc ();
    gsl_spline *spline2 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline2, k_t, P2, th_dim);
    gsl_interp_accel *acc3 = gsl_interp_accel_alloc ();
    gsl_spline *spline3 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline3, k_t, P3, th_dim);
    gsl_interp_accel *acc4 = gsl_interp_accel_alloc ();
    gsl_spline *spline4 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline4, k_t, P4, th_dim);
    gsl_interp_accel *acc5 = gsl_interp_accel_alloc ();
    gsl_spline *spline5 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline5, k_t, P5, th_dim);
    gsl_interp_accel *acc6 = gsl_interp_accel_alloc ();
    gsl_spline *spline6 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline6, k_t, P6, th_dim);
    gsl_interp_accel *acc7 = gsl_interp_accel_alloc ();
    gsl_spline *spline7 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline7, k_t, P7, th_dim);
    gsl_interp_accel *acc8 = gsl_interp_accel_alloc ();
    gsl_spline *spline8 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline8, k_t, P8, th_dim);
    gsl_interp_accel *acc9 = gsl_interp_accel_alloc ();
    gsl_spline *spline9 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline9, k_t, P9, th_dim);
    gsl_interp_accel *acc10 = gsl_interp_accel_alloc ();
    gsl_spline *spline10 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline10, k_t, P10, th_dim);
    gsl_interp_accel *acc11 = gsl_interp_accel_alloc ();
    gsl_spline *spline11 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline11, k_t, P11, th_dim);
    gsl_interp_accel *acc12 = gsl_interp_accel_alloc ();
    gsl_spline *spline12 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline12, k_t, P12, th_dim);
    gsl_interp_accel *acc13 = gsl_interp_accel_alloc ();
    gsl_spline *spline13 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline13, k_t, P13, th_dim);
    gsl_interp_accel *acc14 = gsl_interp_accel_alloc ();
    gsl_spline *spline14 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline14, k_t, P14, th_dim);
    gsl_interp_accel *acc15 = gsl_interp_accel_alloc ();
    gsl_spline *spline15 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline15, k_t, P15, th_dim);
    gsl_interp_accel *acc16 = gsl_interp_accel_alloc ();
    gsl_spline *spline16 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline16, k_t, P16, th_dim);
    gsl_interp_accel *acc17 = gsl_interp_accel_alloc ();
    gsl_spline *spline17 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline17, k_t, P17, th_dim);
    gsl_interp_accel *acc18 = gsl_interp_accel_alloc ();
    gsl_spline *spline18 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline18, k_t, P18, th_dim);
    gsl_interp_accel *acc19 = gsl_interp_accel_alloc ();
    gsl_spline *spline19 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline19, k_t, P19, th_dim);
    gsl_interp_accel *acc20 = gsl_interp_accel_alloc ();
    gsl_spline *spline20 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline20, k_t, P20, th_dim);
    gsl_interp_accel *acc21 = gsl_interp_accel_alloc ();
    gsl_spline *spline21 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline21, k_t, P21, th_dim);
    gsl_interp_accel *acc22 = gsl_interp_accel_alloc ();
    gsl_spline *spline22 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline22, k_t, P22, th_dim);
    gsl_interp_accel *acc23 = gsl_interp_accel_alloc ();
    gsl_spline *spline23 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline23, k_t, P23, th_dim);
    gsl_interp_accel *acc24 = gsl_interp_accel_alloc ();
    gsl_spline *spline24 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline24, k_t, P24, th_dim);
    gsl_interp_accel *acc25 = gsl_interp_accel_alloc ();
    gsl_spline *spline25 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline25, k_t, P25, th_dim);
    gsl_interp_accel *acc26 = gsl_interp_accel_alloc ();
    gsl_spline *spline26 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline26, k_t, P26, th_dim);
    gsl_interp_accel *acc27 = gsl_interp_accel_alloc ();
    gsl_spline *spline27 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline27, k_t, P27, th_dim);
    gsl_interp_accel *acc28 = gsl_interp_accel_alloc ();
    gsl_spline *spline28 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline28, k_t, P28, th_dim);
    gsl_interp_accel *acc29 = gsl_interp_accel_alloc ();
    gsl_spline *spline29 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline29, k_t, P29, th_dim);
    gsl_interp_accel *acc30 = gsl_interp_accel_alloc ();
    gsl_spline *spline30 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline30, k_t, P30, th_dim);
    gsl_interp_accel *acc31 = gsl_interp_accel_alloc ();
    gsl_spline *spline31 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline31, k_t, P31, th_dim);
    gsl_interp_accel *acc32 = gsl_interp_accel_alloc ();
    gsl_spline *spline32 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline32, k_t, P32, th_dim);
    gsl_interp_accel *acc33 = gsl_interp_accel_alloc ();
    gsl_spline *spline33 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline33, k_t, P33, th_dim);
    gsl_interp_accel *acc34 = gsl_interp_accel_alloc ();
    gsl_spline *spline34 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline34, k_t, P34, th_dim);
    gsl_interp_accel *acc35 = gsl_interp_accel_alloc ();
    gsl_spline *spline35 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline35, k_t, P35, th_dim);
    gsl_interp_accel *acc36 = gsl_interp_accel_alloc ();
    gsl_spline *spline36 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline36, k_t, P36, th_dim);
    gsl_interp_accel *acc37 = gsl_interp_accel_alloc ();
    gsl_spline *spline37 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline37, k_t, P37, th_dim);
    gsl_interp_accel *acc38 = gsl_interp_accel_alloc ();
    gsl_spline *spline38 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline38, k_t, P38, th_dim);
    gsl_interp_accel *acc39 = gsl_interp_accel_alloc ();
    gsl_spline *spline39 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline39, k_t, P39, th_dim);
    gsl_interp_accel *acc40 = gsl_interp_accel_alloc ();
    gsl_spline *spline40 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    gsl_spline_init (spline40, k_t, P40, th_dim);
    //gsl_interp_accel *acc41 = gsl_interp_accel_alloc ();
    //gsl_spline *spline41 = gsl_spline_alloc (gsl_interp_cspline, th_dim);
    //gsl_spline_init (spline41, k_t, P41, th_dim);

    
    
    for (i=0; i<kknl_dim; i++)
    {
        struct f_params3 f_par ={alpa,alpe,s8_sc,b1,b2,F,sP,kknl[i],multip,acc1,spline1,acc2,spline2,acc3,spline3,acc4,spline4,acc5,spline5,acc6,spline6,acc7,spline7,acc8,spline8,acc9,spline9,acc10,spline10,acc11,spline11,acc12,spline12,acc13,spline13,acc14,spline14,acc15,spline15,acc16,spline16,acc17,spline17,acc18,spline18,acc19,spline19,acc20,spline20,acc21,spline21,acc22,spline22,acc23,spline23,acc24,spline24,acc25,spline25,acc26,spline26,acc27,spline27,acc28,spline28,acc29,spline29,acc30,spline30,acc31,spline31,acc32,spline32,acc33,spline33,acc34,spline34,acc35,spline35,acc36,spline36,acc37,spline37,acc38,spline38,acc39,spline39,acc40,spline40};
        
        hcubature(1,P024_int,&f_par,1,xmin,xmax,0, 0, 5e-6, 5e-6, &val, &err);
        Pout[i] = val;
    }    

}
