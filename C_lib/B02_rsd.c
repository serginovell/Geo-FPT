#include "B02_rsd.h"
#include "stdio.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include "cubature.h"

/*
 * GEO-FPT
 * All rights reserved
 * Author: Sergi Novell Masot
 * Date: 28th March 2023
 * email: sergi.novell@icc.ub.edu
 * 
*/


//Cosine ka-kb
double cosab(double ka, double kb, double kc)
{
    return ((kc*kc - ka*ka - kb*kb)/(2*ka*kb));
}


//SPT F2 and G2 kernels
double f2_ker(double ka, double kb, double kc)
{
    double cab = cosab(ka,kb,kc);
    return 5./7. +0.5*cab*(ka/kb + kb/ka) + 2./7.*pow(cab,2.0);
    
}
//interpolate kernels across scale factor a(z)
double* interpol_ker(double a, double f1[], double f2[], double f3[]) {
  static double f[5];
  int i;
  if (a>=0.5) { //if the redshift is smaller or equal than 1, obtain the line between z=0.5 and z=1
    for (i = 0;i<5;i++){
      f[i] = f2[i]+(f1[i]-f2[i])*(a-0.5)/(2.0/3.0-0.5);
    }
  }else{ //if redshift is bigger than 1, obtain the line between z=1 and 2
    for (i = 0;i<5;i++){
      f[i] = f2[i]+(f3[i]-f2[i])*(a-0.5)/(1.0/3.0-0.5);
    }
  }
  return f;
}

double g2_ker(double ka, double kb, double kc)
{
    double cab = cosab(ka,kb,kc);

    return 3./7. +0.5*cab*(ka/kb + kb/ka) + 4./7.*pow(cab,2.0);}

// Redshift-space SPT kernels Z1 and Z2
double z1_ker(double mu, double *cosm_par)
{
    double b1 = cosm_par[4], ff = cosm_par[1];
    return b1 + ff * pow(mu,2.);
}

double z2_ker(double ka, double kb, double kc,double fkern, double gkern, double mua, double mub, double *cosm_par)
{
    double cab = cosab(ka,kb,kc);
    double b1 = cosm_par[4], ff = cosm_par[1], b2 =cosm_par[5];

    double ksq = sqrt(ka*ka + kb*kb + 2.0*ka*kb*cab); //modulus of vector sum k1+k2
    double mu12 = (ka*mua + kb*mub)/ksq;
    double bs = -4.0/7.0*(b1 - 1.0);
    double s2 = pow(cab,2) - 1.0/3.0; //S_2 kernel

    double b1_terms = b1 * (fkern +  (0.5*ff*mu12*ksq * (mua/ka + mub/kb)));
    double g_term=0.,fsq_term=0.,b_terms=0.;

    g_term = ff*pow(mu12,2.0)*gkern;
    fsq_term = (0.5 * ff*ff * mu12 * ksq * mua * mub * (mub/ka + mua/kb));
    b_terms = 0.5*(b2 + bs*s2);

    return b1_terms + g_term + fsq_term + b_terms ;
}
//GEO-FPT factor multiplying Z2_SPT to obtain Z2_GEO
double geo_fac(double ka, double kb, double kc, double *af,double hh)
{
    double cosmed=0.,cosmin=0.,cosmax=0.,kmax=0.,kmed=0.,kmin=0.;
    double perim  = (ka + kb + kc)/2., area = sqrt(perim*(perim-ka)*(perim-kb)*(perim-kc))/(hh*hh*0.001);
    if (ka>=kb && ka >= kc){
        kmax=ka;
        if (kb>=kc){kmed=kb; kmin=kc;}
        else       {kmed=kc; kmin=kb;}}
    if (kb>=ka && kb >= kc){
        kmax=kb;
        if (ka>=kc){kmed=ka; kmin=kc;}
        else       {kmed=kc; kmin=ka;}}
    if (kc>=ka && kc >= kb){
        kmax=kc;
        if (ka>=kb){kmed=ka; kmin=kb;}
        else       {kmed=kb; kmin=ka;}}

    cosmax = (pow(kmed,2)+pow(kmin,2)-pow(kmax,2))/(2.*kmed*kmin);
    cosmed = (pow(kmax,2)+pow(kmin,2)-pow(kmed,2))/(2.*kmax*kmin);
    cosmin = (pow(kmax,2)+pow(kmed,2)-pow(kmin,2))/(2.*kmax*kmed);

    double extra = af[0]+af[1]*cosmed/cosmin+af[2]*cosmax/cosmin +af[3]*area+af[4]*area*area;

    return extra;
}


//Integrand
int bkeff_r(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval )
{

    struct input * fp = (struct input *)fdata;
    //load parameters
    double *tr  = fp->tr;
    double *cosm_par = fp->cosm_par;
    double *pk_in = fp->pk_in;
    double sig_fog = fp->sig_fog;
    gsl_interp_accel *acc_me = fp->acc_me;
    gsl_spline *spline_me    = fp->spline_me;
    double *af     = fp->af;
    int mp  = fp->mp;

    //remove warning
    (void)ndim;
    (void)fdim;

    double alpa = cosm_par[2];
    double alpe = cosm_par[3];
    double Fsq = 1./pow(alpa/alpe,2.);

    double ka_m = tr[0];
    double kb_m = tr[1];
    double kc_m = tr[2];
    double pka = pk_in[0];
    double pkb = pk_in[1];
    double pkc = pk_in[2];

    double cab_m = cosab(ka_m,kb_m,kc_m);
    //Integration angles
    double mua_m = x[0];
    double phi = x[1];
    
    //obtain the cosines respect to los of kb and kc
    double mub_m = mua_m*cab_m - sqrt((1.0-pow(mua_m,2.0))*(1.0-pow(cab_m,2.0)))*cos(phi);
    double muc_m = (-ka_m*mua_m-kb_m*mub_m)/kc_m;
    //Apply AP transformation for k's
    double ka  = ka_m * sqrt(1.+pow(mua_m,2.)*(Fsq-1.))/alpe;
    double kb  = kb_m * sqrt(1.+pow(mub_m,2.)*(Fsq-1.))/alpe;
    double kc  = kc_m * sqrt(1.+pow(muc_m,2.)*(Fsq-1.))/alpe;
    double hh = 1.; //if working in Mpc, set hh=0.6711 (The h of Quijote fiducial simulations)

    if(kb+ka-kc<hh*1.1*2*M_PI/1000.||ka+kc-kb<hh*1.1*2*M_PI/1000.||(kb+kc)-ka<hh*1.1*2*M_PI/1000.){
        //discard if the new triangle has a closure condition smaller than Delta k. If Delta k is higher this can be changed, although the difference should be minimal
	fval[0]=0;
	}else{
    //Apply AP transformation for mu's
    double mua = mua_m*alpe/(alpa*sqrt(1.+pow(mua_m,2.)*(Fsq-1.)) );
    double mub = mub_m*alpe/(alpa*sqrt(1.+pow(mub_m,2.)*(Fsq-1.)) );
    double muc = muc_m*alpe/(alpa*sqrt(1.+pow(muc_m,2.)*(Fsq-1.)) );
    
    double eff_fact = geo_fac(ka,kb,kc,af,hh);
    double D_fog= 1./pow(1.0 + 0.5 * pow(pow(ka*mua, 2.0)+pow(kb*mub, 2.0)+pow(kc*muc, 2.0),2.0) * pow(sig_fog/hh,4.0), 2.0); //FoG damping factor

    //calculate kernels
    double z1_1 = z1_ker(mua,cosm_par);
    double z1_2 = z1_ker(mub,cosm_par);
    double z1_3 = z1_ker(muc,cosm_par);


    double f2k_12 = f2_ker(ka,kb,kc);
    double f2k_23 = f2_ker(kc,kb,ka);
    double f2k_13 = f2_ker(ka,kc,kb);

    double g2k_12 = g2_ker(ka,kb,kc);
    double g2k_23 = g2_ker(kc,kb,ka);
    double g2k_13 = g2_ker(ka,kc,kb);

    double z2_12 = z2_ker(ka,kb,kc,f2k_12,g2k_12,mua,mub,cosm_par)*eff_fact;
    double z2_23 = z2_ker(kb,kc,ka,f2k_23,g2k_23,mub,muc,cosm_par)*eff_fact;
    double z2_13 = z2_ker(ka,kc,kb,f2k_13,g2k_13,mua,muc,cosm_par)*eff_fact;
    
    //Legendre polynomial, depending on the multipole (mp)
    double leg = 1;
    if (mp == 1){leg = 5*(3.*pow(mua,2)-1.)/2.;}
    if (mp == 2){leg = 5*(3.*pow(mub,2)-1.)/2.;}
    if (mp == 3){leg = 5*(3.*pow(muc,2)-1.)/2.;}

        
    //Interpolate power spectrum
    pka = pow(10.0, gsl_spline_eval(spline_me, log10(ka), acc_me));
    pkb = pow(10.0, gsl_spline_eval(spline_me, log10(kb), acc_me));
    pkc = pow(10.0, gsl_spline_eval(spline_me, log10(kc), acc_me));

    fval[0]  = leg*D_fog*(z1_1*z1_2*z2_12*pka*pkb + z1_3*z1_2*z2_23*pkc*pkb + z1_1*z1_3*z2_13*pka*pkc)/(2.*M_PI*pow(alpa,2.)*pow(alpe,4.));
    }
    return 0;
}

//Integration
extern void ext_bk_mp(double **tr,double **tr2, double **tr3, double **tr4, double *log_km, double *log_pkm, double *cosm_par, double redshift, int fit_full,int kp_dim,  
int num_tr, int num_tr2,int num_tr3,int num_tr4, double *bk_mipj_ar)
{
    
    //interpolate the calibrated f1,...,f5 parameters at the redshift
    double a_t = 1./(1.+redshift);
    static double* af;
    if(fit_full==1){
        double f1[5]={1.00334041e+00, -3.97587153e-03,  2.39676798e-02, -5.68098245e-02,  1.32544359e-02};
        double f2[5]={1.01798746e+00, -4.10578402e-03,  1.48808466e-02, -5.46660187e-02,1.11441504e-02};
        double f3[5]={1.03745852e+00,  1.95019784e-03,  5.56349564e-03, -4.83637541e-02, 8.06689225e-03};
        af = interpol_ker(a_t, f1,f2, f3);
    }else{
        double f1[5]={ 1.00321199,  0.00778061,  0.04035794, -0.07043621,  0.01390807};
        double f2[5]={ 1.15615701,  0.01335707, -0.00582238, -0.08951067,  0.01584383};
        double f3[5]={ 1.29787528,  0.01154473, -0.06984416, -0.09486691,  0.01642198};
        af = interpol_ker(a_t, f1,f2, f3);
    }
    
    //printf("%.8f %.8f %.8f %.8f %.8f\n",af[0],af[1],af[2],af[3],af[4]);
    
    //Interpolate Non-linear Pk
    gsl_interp_accel *acc_me = gsl_interp_accel_alloc ();
    gsl_spline *spline_me = gsl_spline_alloc (gsl_interp_cspline, kp_dim);
    gsl_spline_init (spline_me, log_km, log_pkm, kp_dim); 

    double val, err;
    double sig_fog = cosm_par[9];
    double xmin[2] = {-1.0 , 0.0};
    double xmax[2] = {1.0 , 2.0*M_PI};
    int i;
    //monopole
    for (i = 0 ; i < num_tr; i ++)
    {
        double ka= tr[i][0];
        double kb= tr[i][1];
        double kc= tr[i][2];
        //interpolate non-linear pk
        double pka = pow(10.0, gsl_spline_eval(spline_me, log10(ka), acc_me));
        double pkb = pow(10.0, gsl_spline_eval(spline_me, log10(kb), acc_me));
        double pkc = pow(10.0, gsl_spline_eval(spline_me, log10(kc), acc_me));

        double pk_in[3] = {pka,pkb,pkc};
        struct input params = {tr[i], cosm_par, pk_in,sig_fog,acc_me,spline_me,af,0};
        hcubature(1,bkeff_r, &params, 2, xmin, xmax, 0, 0, 2e-4, ERROR_INDIVIDUAL, &val, &err);

        bk_mipj_ar[i] = val;
    }
    
    //B200
    for (i = 0 ; i < num_tr2; i ++)
    {
        double ka= tr2[i][0];
        double kb= tr2[i][1];
        double kc= tr2[i][2];

        //interpolate non-linear pk
        double pka = pow(10.0, gsl_spline_eval(spline_me, log10(ka), acc_me));
        double pkb = pow(10.0, gsl_spline_eval(spline_me, log10(kb), acc_me));
        double pkc = pow(10.0, gsl_spline_eval(spline_me, log10(kc), acc_me));
        
        double pk_in[3] = {pka,pkb,pkc};
         
        struct input params = {tr2[i], cosm_par, pk_in,sig_fog,acc_me,spline_me,af,1};
        hcubature(1,bkeff_r, &params, 2, xmin, xmax, 0, 0, 5e-4, ERROR_INDIVIDUAL, &val, &err);
            
        bk_mipj_ar[i+num_tr] = val;
    }
    //B020
    for (i = 0 ; i < num_tr3; i ++)
    {
        double ka= tr3[i][0];
        double kb= tr3[i][1];
        double kc= tr3[i][2];


        //interpolate non-linear pk  
        double pka = pow(10.0, gsl_spline_eval(spline_me, log10(ka), acc_me));
        double pkb = pow(10.0, gsl_spline_eval(spline_me, log10(kb), acc_me));
        double pkc = pow(10.0, gsl_spline_eval(spline_me, log10(kc), acc_me));
        
        double pk_in[3] = {pka,pkb,pkc};
         
        struct input params = {tr3[i], cosm_par, pk_in,sig_fog,acc_me,spline_me,af,2};
        hcubature(1,bkeff_r, &params, 2, xmin, xmax, 0, 0, 5e-4, ERROR_INDIVIDUAL, &val, &err);
            
        bk_mipj_ar[i+num_tr+num_tr2] = val;
    }

    //B002
    for (i = 0 ; i < num_tr4; i ++)
    {
        double ka= tr4[i][0];
        double kb= tr4[i][1];
        double kc= tr4[i][2];

        //interpolate non-linear pk    
        double pka = pow(10.0, gsl_spline_eval(spline_me, log10(ka), acc_me));
        double pkb = pow(10.0, gsl_spline_eval(spline_me, log10(kb), acc_me));
        double pkc = pow(10.0, gsl_spline_eval(spline_me, log10(kc), acc_me));
        
        double pk_in[3] = {pka,pkb,pkc};
         
        struct input params = {tr4[i], cosm_par,pk_in,sig_fog,acc_me,spline_me,af,3};
        hcubature(1,bkeff_r, &params, 2, xmin, xmax, 0, 0, 5e-4, ERROR_INDIVIDUAL, &val, &err);
            
        bk_mipj_ar[i+num_tr+num_tr2+num_tr3] = val;
     }

    gsl_spline_free (spline_me);
    gsl_interp_accel_free (acc_me);

}
