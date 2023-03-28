#include "r_space.h"
#include "stdio.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include "cubature.h"

//Cosine ka-kb
double cosab(double ka, double kb, double kc)
{
    return ((kc*kc - ka*ka - kb*kb)/(2*ka*kb));
}

//Fitting functions
double a_eff(double nn, double qq, double s8, double *aa)
{
    double q3  = (4.0 - pow(2.0, nn))/(1.0 + pow(2.0, nn + 1.0));
    double num = 1.0 + pow(s8, aa[5]) * pow(0.7 * q3, 0.5) * pow(qq * aa[0], nn + aa[1]);
    double den = 1.0 + pow(qq * aa[0], nn + aa[1]);

    return num/den;
}
double b_eff(double nn, double qq, double * aa)
{
    double num = 1.0 + 0.2 * aa[2] * (nn + 3.0) * pow(qq * aa[6], nn + 3.0 + aa[7]);
    double den = 1.0 + pow(qq * aa[6], nn + 3.5 + aa[7]);
    return num/den;
}

double c_eff(double nn, double qq, double * aa)
{
    double num = 1.0 + 4.5 * aa[3] / (1.5 + pow(nn + 3.0, 4.0)) * pow(qq * aa[4], nn + 3.0 + aa[8]);
    double den = 1.0 + pow(qq * aa[4], nn + 3.5 + aa[8]);

    return num/den;
}

double b_eff_n(double nn, double qq, double * aa)
{
    double num = 1.0 + 0.2 * aa[0] * (nn + 3.0) * pow(qq * aa[3], nn + 3.0 + aa[4]);
    double den = 1.0 + pow(qq * aa[3], nn + 3.5 + aa[4]);
    return num/den;
}

double c_eff_n(double nn, double qq, double * aa)
{
    double num = 1.0 + 4.5 * aa[1] / (1.5 + pow(nn + 3.0, 4.0)) * pow(qq * aa[2], nn + 3.0 + aa[5]);
    double den = 1.0 + pow(qq * aa[2], nn + 3.5 + aa[5]);

    return num/den;
}
double ome(double z, double omm)
{
    double Om0 = omm;
    double om  = Om0*pow(1.+z,3)/(Om0*pow(1.+z,3)+1.-Om0);
    double t   = 1.0/om - 1.0;
    double L   = log(pow(1+t,0.5)-pow(t,0.5));
    double D   = 1. +3./t + 3*pow(1.+t,0.5)/pow(t,1.5)*L;
    double sq  = pow((1.+t)/t,0.5);
    double aux = 1 + sq*L + 0.5*pow(sq+L/t,2);
    return 0.5-1./(4.*D*D)-9./(4.*t*D*D)*aux;
    }
//Effective F and G kernels
double f2_ker(double ka, double kb, double kc, double na, double nb, double knl,  double s8, double mod, double *af)
{
    double cab = cosab(ka,kb,kc);
    double out = 0.;
    //SPT
    if(mod==3){out = 5./7. +0.5*cab*(ka/kb + kb/ka) + 2./7.*pow(cab,2.0);}
    //EFF
    if(mod==1){
        double ai = a_eff(na, ka/knl, s8, af);
        double aj = a_eff(nb, kb/knl, s8, af);
        double bi = b_eff(na, ka/knl, af);
        double bj = b_eff(nb, kb/knl, af);
        double ci = c_eff(na, ka/knl, af);
        double cj = c_eff(nb, kb/knl, af);

        out = 5./7.*ai*aj + 0.5*cab*bi*bj*(ka/kb + kb/ka) + 2./7.*pow(cab,2.0)*ci*cj;
        }


    return out;
}

double g2_ker(double ka, double kb, double kc, double na, double nb, double knl,  double s8, double mod, double *ag)
{
    double cab = cosab(ka,kb,kc);
    double out = 0.;
    //SPT
    if(mod==3){out = 3./7. +0.5*cab*(ka/kb + kb/ka) + 4./7.*pow(cab,2.0);}
    //EFF
    if(mod==1){
        double ai = a_eff(na, ka/knl, s8, ag);
        double aj = a_eff(nb, kb/knl, s8, ag);
        double bi = b_eff(na, ka/knl, ag);
        double bj = b_eff(nb, kb/knl, ag);
        double ci = c_eff(na, ka/knl, ag);
        double cj = c_eff(nb, kb/knl, ag);
        out = 3./7.*ai*aj + 0.5*cab*bi*bj*(ka/kb + kb/ka) + 4./7.*pow(cab,2.0)*ci*cj;}

    return out;
}


// Redshift-space kernels
double z1_ker(double mu, double *cosm_par)
{
    double b1 = cosm_par[1], ff = cosm_par[2];
    return b1 + ff * pow(mu,2.);
}

double z2_ker(double ka, double kb, double kc,double fkern, double gkern, double mua, double mub, double *cosm_par,int mat)
{
    double cab = cosab(ka,kb,kc);
    double b1 = cosm_par[1], ff = cosm_par[2], b2 =cosm_par[3];

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

double geo_fac(double ka, double kb, double kc, double *af,double *ag, int ell,double hh)
{
    double cosmed=0.,cosmin=0.,cosmax=0.,kmax=0.,kmed=0.,kmin=0.;
    double perim  = (ka + kb + kc)/2., sqarea = sqrt(perim*(perim-ka)*(perim-kb)*(perim-kc))/(hh*hh*0.001);
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

//    double extra = af[0]+ell*ag[0]+(af[1]+ell*ag[1])*cosmed/cosmin+(af[2]+ell*ag[2])*pow(cosmed/cosmin,2)+(af[3]+ell*ag[3])*cosmax/cosmin+(af[4]+ell*ag[4])*pow(cosmax/cosmin,2) +(af[5]+ell*ag[5])*(sqarea)+(af[6]+ell*ag[6])*sqarea*sqarea;
    double extra = af[0]+ell*ag[0]+(af[1]+ell*ag[1])*cosmed/cosmin+(af[2]+ell*ag[2])*cosmax/cosmin +(af[3]+ell*ag[3])*(sqarea)+(af[4]+ell*ag[4])*sqarea*sqarea;

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
    gsl_interp_accel *acc_n = fp->acc_n;
    gsl_spline *spline_n    = fp->spline_n;
    gsl_interp_accel *acc_me = fp->acc_me;
    gsl_spline *spline_me    = fp->spline_me;
    double knl              = fp->knl;
    double D_z     = fp->D_z;
    int mod        = fp->mod;
    double *af     = fp->af;
    double *ag     = fp->ag;
    int mat = fp->mat;
    int mp  = fp->mp;

    //remove warning
    (void)ndim;
    (void)fdim;

    double s8 = cosm_par[0];
    double alpa = cosm_par[5];
    double alpe = cosm_par[6];
    double Fsq = 1./pow(alpa/alpe,2.);

    double ka_m = tr[0];
    double kb_m = tr[1];
    double kc_m = tr[2];
    double pka = pk_in[0];
    double pkb = pk_in[1];
    double pkc = pk_in[2];
    double na=0.; 
    double nb=0.;
    double nc=0.;
    double cab_m = cosab(ka_m,kb_m,kc_m);
    double mua_m = x[0];
    double phi = x[1];
    double mub_m = mua_m*cab_m - sqrt((1.0-pow(mua_m,2.0))*(1.0-pow(cab_m,2.0)))*cos(phi);
    double muc_m = (-ka_m*mua_m-kb_m*mub_m)/kc_m;

    double ka  = ka_m * sqrt(1.+pow(mua_m,2.)*(Fsq-1.))/alpe;
    double kb  = kb_m * sqrt(1.+pow(mub_m,2.)*(Fsq-1.))/alpe;
    double kc  = kc_m * sqrt(1.+pow(muc_m,2.)*(Fsq-1.))/alpe;
    double hh = 1.;

    if(kb+ka-kc<hh*1.1*2*M_PI/1000.||ka+kc-kb<hh*1.1*2*M_PI/1000.||(kb+kc)-ka<hh*1.1*2*M_PI/1000.){
	fval[0]=0;
	}else{
    double mua = mua_m*alpe/(alpa*sqrt(1.+pow(mua_m,2.)*(Fsq-1.)) );
    double mub = mub_m*alpe/(alpa*sqrt(1.+pow(mub_m,2.)*(Fsq-1.)) );
    double muc = muc_m*alpe/(alpa*sqrt(1.+pow(muc_m,2.)*(Fsq-1.)) );
    
    double eff_fact = 1.;
    double ell =0.;

    if(mod==3){eff_fact=geo_fac(ka,kb,kc,af,ag,ell,hh);}
    double D_fog;
    // Fingers-of-God damping factor


    D_fog = 1./pow(1.0 + 0.5 * pow(pow(ka*mua, 2.0)+pow(kb*mub, 2.0)+pow(kc*muc, 2.0),2.0) * pow(sig_fog/hh,4.0), 2.0);

    //calculate kernels
    double z1_1 = z1_ker(mua,cosm_par);
    double z1_2 = z1_ker(mub,cosm_par);
    double z1_3 = z1_ker(muc,cosm_par);
    double z2_12, z2_13,z2_23;
    
    
    pka = pow(10.0, gsl_spline_eval(spline_me, log10(ka), acc_me));
    pkb = pow(10.0, gsl_spline_eval(spline_me, log10(kb), acc_me));
    pkc = pow(10.0, gsl_spline_eval(spline_me, log10(kc), acc_me));
    if(mod==1){
    na = gsl_spline_eval(spline_n, log10(ka), acc_n);
    nb = gsl_spline_eval(spline_n, log10(kb), acc_n);
    nc = gsl_spline_eval(spline_n, log10(kc), acc_n);
    }
    double f2k_12 = f2_ker(ka,kb,kc,na,nb,knl,s8,mod,af);
    double f2k_23 = f2_ker(kc,kb,ka,nc,nb,knl,s8,mod,af);
    double f2k_13 = f2_ker(ka,kc,kb,na,nc,knl,s8,mod,af);

    double g2k_12 = g2_ker(ka,kb,kc,na,nb,knl,s8,mod,ag);
    double g2k_23 = g2_ker(kc,kb,ka,nc,nb,knl,s8,mod,ag);
    double g2k_13 = g2_ker(ka,kc,kb,na,nc,knl,s8,mod,ag);

    z2_12 = z2_ker(ka,kb,kc,f2k_12,g2k_12,mua,mub,cosm_par,mat)*eff_fact;
    z2_23 = z2_ker(kb,kc,ka,f2k_23,g2k_23,mub,muc,cosm_par,mat)*eff_fact;
    z2_13 = z2_ker(ka,kc,kb,f2k_13,g2k_13,mua,muc,cosm_par,mat)*eff_fact;
    
    //Legendre polynomial, depending on the multipole (mp)
    double leg = 1;
    if (mp == 1){leg = 5*(3.*pow(mua,2)-1.)/2.;}
    if (mp == 2){leg = 5*(3.*pow(mub,2)-1.)/2.;}
    if (mp == 3){leg = 5*(3.*pow(muc,2)-1.)/2.;}


    if(mat==0){fval[0]  = leg*D_fog*(z1_1*z1_2*z2_12*pka*pkb + z1_3*z1_2*z2_23*pkc*pkb + z1_1*z1_3*z2_13*pka*pkc)/(2.*M_PI*pow(alpa,2.)*pow(alpe,4.));}

    if(mat==1){fval[0] = (z1_1*z1_2*z2_12*pka*pkb + z1_3*z1_2*z2_23*pkc*pkb + z1_1*z1_3*z2_13*pka*pkc)/(2.*M_PI*pow(alpa,2.)*pow(alpe,4.));} //real space

    }
    return 0;
}

//Integration
extern void ext_bk_mp(double **tr,double **tr2, double **tr3, double **tr4, double *log_km, double *log_pkm, double *log_kh_n, double *n_smooth, double *cosm_par, double *af, double *ag, double *af_mup, double *ag_mup,int tr_dim,  
int kh_n_dim, int num_tr, int num_tr2,int num_tr3,int num_tr4, double *bk_mipj_ar, double knl, double D_z, int mod,double sig_fog,double sig_fog2, int multip, int mat)
{
    double s8 = cosm_par[0];
    //Interpolate smoothed slope of linear Pk
    gsl_interp_accel *acc_n = gsl_interp_accel_alloc ();
    gsl_spline *spline_n = gsl_spline_alloc (gsl_interp_cspline, kh_n_dim);
    gsl_spline_init (spline_n, log_kh_n, n_smooth, kh_n_dim);
    //Interpolate measured Pk
    gsl_interp_accel *acc_me = gsl_interp_accel_alloc ();
    gsl_spline *spline_me = gsl_spline_alloc (gsl_interp_cspline, tr_dim);
    gsl_spline_init (spline_me, log_km, log_pkm, tr_dim); 

    double val, err;
    double xmin[2] = {-1.0 , 0.0};
    double xmax[2] = {1.0 , 2.0*M_PI};
    int i;
    //monopole
    for (i = 0 ; i < num_tr; i ++)
    {
        double ka= tr[i][0];
        double kb= tr[i][1];
        double kc= tr[i][2];
        //interpolate linear or measured pk and smoothed n    
        double pka = pow(10.0, gsl_spline_eval(spline_me, log10(ka), acc_me));
        double pkb = pow(10.0, gsl_spline_eval(spline_me, log10(kb), acc_me));
        double pkc = pow(10.0, gsl_spline_eval(spline_me, log10(kc), acc_me));

        double pk_in[3] = {pka,pkb,pkc};
        struct input params = {tr[i], cosm_par, pk_in,sig_fog,1.,acc_n,spline_n,acc_n,spline_n,acc_me,spline_me,knl,D_z,mod,af,ag,mat,0};
        hcubature(1,bkeff_r, &params, 2, xmin, xmax, 0, 0, 2e-4, ERROR_INDIVIDUAL, &val, &err);

        bk_mipj_ar[i] = val;
    }
    
    if(multip>2.5){
        //B200
        for (i = 0 ; i < num_tr2; i ++)
        {
            double ka= tr2[i][0];
            double kb= tr2[i][1];
            double kc= tr2[i][2];

            //interpolate linear or measured pk and smoothed n    
            double pka = pow(10.0, gsl_spline_eval(spline_me, log10(ka), acc_me));
            double pkb = pow(10.0, gsl_spline_eval(spline_me, log10(kb), acc_me));
            double pkc = pow(10.0, gsl_spline_eval(spline_me, log10(kc), acc_me));
        
            double pk_in[3] = {pka,pkb,pkc};
         
            struct input params = {tr2[i], cosm_par, pk_in,sig_fog,1.,acc_n,spline_n,acc_n,spline_n,acc_me,spline_me,knl,D_z,mod,af,ag,mat,1};
            hcubature(1,bkeff_r, &params, 2, xmin, xmax, 0, 0, 5e-4, ERROR_INDIVIDUAL, &val, &err);
            
            bk_mipj_ar[i+num_tr] = val;
        }
        //B020
        if(multip>3.5){
        for (i = 0 ; i < num_tr3; i ++)
        {
            double ka= tr3[i][0];
            double kb= tr3[i][1];
            double kc= tr3[i][2];


            //interpolate linear or measured pk and smoothed n    
            double pka = pow(10.0, gsl_spline_eval(spline_me, log10(ka), acc_me));
            double pkb = pow(10.0, gsl_spline_eval(spline_me, log10(kb), acc_me));
            double pkc = pow(10.0, gsl_spline_eval(spline_me, log10(kc), acc_me));
        
            double pk_in[3] = {pka,pkb,pkc};
         
            struct input params = {tr3[i], cosm_par, pk_in,sig_fog,1.,acc_n,spline_n,acc_n,spline_n,acc_me,spline_me,knl,D_z,mod,af,ag,mat,2};
            hcubature(1,bkeff_r, &params, 2, xmin, xmax, 0, 0, 5e-4, ERROR_INDIVIDUAL, &val, &err);
            
            bk_mipj_ar[i+num_tr+num_tr2] = val;
        }
        }
        //B002
        if(multip<-1.1||multip>7.5){
        for (i = 0 ; i < num_tr4; i ++)
        {
            double ka= tr4[i][0];
            double kb= tr4[i][1];
            double kc= tr4[i][2];

            //interpolate linear or measured pk and smoothed n    
            double pka = pow(10.0, gsl_spline_eval(spline_me, log10(ka), acc_me));
            double pkb = pow(10.0, gsl_spline_eval(spline_me, log10(kb), acc_me));
            double pkc = pow(10.0, gsl_spline_eval(spline_me, log10(kc), acc_me));
        
            double pk_in[3] = {pka,pkb,pkc};
         
            struct input params = {tr4[i], cosm_par,pk_in,sig_fog2,1.,acc_n,spline_n,acc_n,spline_n,acc_me,spline_me,knl,D_z,mod,af_mup,ag_mup,mat,3};
            hcubature(1,bkeff_r, &params, 2, xmin, xmax, 0, 0, 5e-3, ERROR_INDIVIDUAL, &val, &err);
            
            bk_mipj_ar[i+num_tr+num_tr2+num_tr3] = val;
        }}
        
        }
      

    gsl_spline_free (spline_n);
    gsl_interp_accel_free (acc_n);
    gsl_spline_free (spline_me);
    gsl_interp_accel_free (acc_me);

}
