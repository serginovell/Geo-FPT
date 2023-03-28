#include "davide.h"
#include "stdio.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include "cubature.h"

double tab_f(double ka,double kb, double kc)
{
    double arg  = (pow(ka,2.0) + pow(kb,2.0) - pow(kc,2.0))/(2.0 * ka * kb);
    double res  = M_PI - acos(arg);

    if (isnan(res)){printf("arg %.5f ka %.5f kb %.5f kc %.5f \n",arg,ka,kb,kc);}

    return res;
}
double s2spt(const double ctab)
{
    return pow(ctab,2.0) - 1.0/3.0;
}

double f2spt(const double ka, const double kb, const double ctab)
{return 5.0 / 7.0 + 0.5 * ctab * (ka/kb + kb/ka) + (2.0 / 7.0) * pow(ctab,2.0);}


double g2spt(const double ka, const double kb, const double ctab)
{return 3.0 / 7.0 + 0.5 * ctab * (ka/kb + kb/ka) + (4.0 / 7.0) * pow(ctab,2.0);}

double z1spt(const double mu, const double *par_ar)
{return par_ar[0]+par_ar[3]*mu*mu;}

double fitz2spt(const double ka, const double kb, const double mua, const double mub,
                const double ctab, const double *par_ar,const double *fit_ar,const double ell)
{
    double b1 = par_ar[0], b2 = par_ar[1], fg = par_ar[3],
            b1_term,b2_term,fg_term,f2_term,res,extra =1.,
            ks = sqrt(pow(ka,2.0) + pow(kb,2.0) + 2.0 * ka * kb * ctab),mu,bs;

    mu = (ka * mua + kb * mub) / ks;

    bs = - 4.0 * (b1 - 1.0)/7.0;

    b1_term = b1 * (f2spt(ka,kb,ctab)
                        + 0.5 * fg * ks * mu * (mub/kb + mua/ka));

    fg_term = fg * pow(mu,2.0) * g2spt(ka,kb,ctab);

    f2_term = 0.5 * pow(fg,2.0) * ks * mu * (mua * mub * (mua/kb + mub/ka));

    b2_term = 0.5 * (b2 + bs * s2spt(ctab));


    double f1 = fit_ar[0], f2 = fit_ar[1], f3 = fit_ar[2], f4 = fit_ar[3], f5 = fit_ar[4], f6 = fit_ar[5], f7 = fit_ar[6],
           g1 = fit_ar[7], g2 = fit_ar[8], g3 = fit_ar[9], g4 = fit_ar[10], g5 = fit_ar[11], g6 = fit_ar[12], g7 = fit_ar[13];
    double kc = sqrt(pow(ka,2)+pow(kb,2)+2.*ka*kb*ctab),
           cosmed=0.,cosmin=0.,kmax=0.,kmed=0.,kmin=0., cosmax=0.,
           perim  = (ka + kb + kc)/2., sqarea = sqrt(perim*(perim-ka)*(perim-kb)*(perim-kc))/0.001;


    if (ka>=kb && ka >= kc)
    {
        kmax=ka;
        if (kb>=kc){kmed=kb; kmin=kc;}
        else       {kmed=kc; kmin=kb;}
    }
    if (kb>=ka && kb >= kc)
    {
        kmax=kb;
        if (ka>=kc){kmed=ka; kmin=kc;}
        else       {kmed=kc; kmin=ka;}
    }
    if (kc>=ka && kc >= kb)
    {
        kmax=kc;
        if (ka>=kb){kmed=ka; kmin=kb;}
        else       {kmed=kb; kmin=ka;}
    }
    cosmax = (pow(kmed,2)+pow(kmin,2)-pow(kmax,2))/(2.*kmed*kmin);
    cosmed = (pow(kmax,2)+pow(kmin,2)-pow(kmed,2))/(2.*kmax*kmin);
    cosmin = (pow(kmax,2)+pow(kmed,2)-pow(kmin,2))/(2.*kmax*kmed);

    extra = (f1 + ell * g1) + (f2 + ell * g2) *      (cosmed/cosmin)
                            + (f3 + ell * g3) *  pow((cosmed/cosmin),2.)
                            + (f4 + ell * g4) *      (cosmax/cosmin)
                            + (f5 + ell * g5) *  pow((cosmax/cosmin),2.)
                            + (f6 + ell * g6) *      sqarea
                            + (f7 + ell * g7) *  pow(sqarea,2.);

    res = extra * (b1_term + fg_term + f2_term + b2_term);
    //printf("ka %.4f kb %.4f kc %.4f cosmax %.4f term1 %.4f term6 %.4f ext %.4f \n", ka,kb,kc, cosmax,f1+ell*g1,(f6+ell*g6)*(sqarea),extra);
    //printf("ka %.4f kb %.4f kc %.4f extra %.4f\n",ka,kb,kc, extra);
    //if(mua<-0.95 && ka<0.05){printf("ka %.4f kb %.4f kc %.4f ext %.4f cosmax %.4f f7 %.4f",ka,kb,kc,extra,cosmax,f7);}
    return res;
}
double new_mufacs(const int idx_mup, const double mua, const double mub, const double muc)
{
    double mup_pref[7];

    mup_pref[0] = 1.;

    mup_pref[1] = 2.5*(3*muc*muc-1.),  mup_pref[2] = muc*muc,  mup_pref[3] = muc*muc;
    //L4
    mup_pref[4] = muc*muc*muc*muc,  mup_pref[5] = muc*muc*muc*muc,   mup_pref[6] = muc*muc*muc*muc;

    return mup_pref[idx_mup];
}

///// Bispectrum
double fog_B( double kpa, double kpb, double kpc, double sb)
{return 1./pow(1. + 0.5 * pow(pow(kpa, 2.) + pow(kpb, 2.) + pow(kpc, 2.),2.) * pow(sb,4.0), 2.);}

int bk_cbrt_gal(unsigned ndim, const double *x, void *params, unsigned fdim, double *fval )
{
    (void)ndim,(void)fdim;
    struct str_bkcbrt * fp = (struct str_bkcbrt*)params;

    //// params ////
    double *kh_ar  = fp->kh_ar, *par_ar = fp->par_ar, *fit_ar = fp->fit_ar,
    alpa = par_ar[4], alpe = par_ar[5], sb = par_ar[6],ff,ffm2,t12,
    mu1,cos_phi,cos_t12,mu2,k1pa,k2pa,k3pa,k1pe,k2pe,k3pe,mu3,
    qq1,qq2,qq3,nu1,nu2,nu3,pq1,pq2,pq3,z1_1,z1_2,z1_3,z2_12,z2_13,z2_23,temp12,temp13,temp23,mup_pr,
    tq12,tq13,tq23,k1,k2,k3,fog_pref=1.,ell=0.,fogk1=1.,fogk2=1.,fogk3=1.;

    int    idx_mup = fp->idx_mup;
    gsl_interp_accel *acc = fp->acc;
    gsl_spline *spline    = fp->spline;

    ff   = alpa/alpe, ffm2 = 1./pow(ff,2.);

    k1 = kh_ar[0], k2 = kh_ar[1], k3 = kh_ar[2];
    //if(mu1<-0.985 && idx_mup==0){printf("ka %.4f kb %.4f kc\n",k1,k2,k3);}
    t12 = tab_f(k1, k2, k3);

    mu1 = x[0], cos_phi = cos(x[1]), cos_t12 = cos(t12);
    mu2 = mu1 * cos_t12 - sqrt(1.0 - pow(mu1,2.0)) * sqrt(1.0 - pow(cos_t12,2.0)) * cos_phi;

    k1pa = k1 * mu1, k2pa = k2 * mu2, k3pa = - k1pa - k2pa;
/*
    k1pe = k1 * sqrt(1 - pow(k1pa/k1,2.0)), k2pe = k2 * sqrt(1 - pow(k2pa/k2,2.0)),k3pe = 0.0;

    if (fabs(k3pa) < k3)
    {k3pe = k3 * sqrt(1 - pow(k3pa/k3,2.0));}

   if (fabs(fabs(k3pa) - k3) < KMINNONZERO)
    {k3pe = 0.;}

    if (fabs(k3pa) > (k3 + KMINNONZERO ) || fabs(k3pe) > (k3 + KMINNONZERO ) ||
        fabs(k2pa) > (k2 + KMINNONZERO ) || fabs(k2pe) > (k2 + KMINNONZERO ) ||
        fabs(k1pa) > (k1 + KMINNONZERO ) || fabs(k1pe) > (k1 + KMINNONZERO ))
    {
        printf("bk check  k1 %.3f, k2 %.3f, k3 %.3f, mu1%.3f, mu2 %.3f, mu3 %.3f, k1pa %.3f, k2pa %.3f, k3pa  %.3f \n",
                k1,k2,k3, mu1, mu2, k3pa/k3, k1pa, k2pa, k3pa);
        fval[0] = 0.0;
    }
    else*/
    
    {
        mu3 = k3pa/k3;

       // compute the fiducial sides and cosines due to AP effect
        qq1 = (k1/alpe) * sqrt(1. + pow(mu1,2.) * (ffm2 - 1.));
        qq2 = (k2/alpe) * sqrt(1. + pow(mu2,2.) * (ffm2 - 1.));
        qq3 = (k3/alpe) * sqrt(1. + pow(mu3,2.) * (ffm2 - 1.));
        nu1 = (mu1/ff)  / sqrt(1. + pow(mu1,2.) * (ffm2 - 1.));
        nu2 = (mu2/ff)  / sqrt(1. + pow(mu2,2.) * (ffm2 - 1.));
        nu3 = (mu3/ff)  / sqrt(1. + pow(mu3,2.) * (ffm2 - 1.));

        //check that the triangle formed by q1, q2, q3 is still closed
        if ((qq1+qq2)>qq3 && (qq2+qq3)>qq1 && (qq1+qq3)>qq2)
        {
             pq1 = gsl_spline_eval(spline, log10(qq1), acc);
             pq2 = gsl_spline_eval(spline, log10(qq2), acc);
             pq3 = gsl_spline_eval(spline, log10(qq3), acc);
             tq12 = tab_f(qq1, qq2, qq3), tq13 = tab_f(qq1, qq3, qq2), tq23 = tab_f(qq2, qq3, qq1);

            if (par_ar[3]>0.00001)
            {
                if (idx_mup>0){ell=1.;}
                fog_pref = fog_B(qq1*nu1,qq2*nu2,qq3*nu3,sb);
                //if ((int)(fit_ar[12])==0){fog_pref = fog_B(qq1*nu1,qq2*nu2,qq3*nu3,fit_ar[10]);}
//                {fogk1 = fog_P_pow(qq1*nu1,fit_ar[10],fit_ar[38],fit_ar[39],ell);
//                 fogk2 = fog_P_pow(qq2*nu2,fit_ar[10],fit_ar[38],fit_ar[39],ell);
//                 fogk3 = fog_P_pow(qq3*nu3,fit_ar[10],fit_ar[38],fit_ar[39],ell);}
//                 {fog_pref = fog_B_atan(qq1*nu1,qq2*nu2,qq3*nu3,fit_ar[10],fit_ar[38],fit_ar[39],ell);}
//                 {fog_pref = fog_B_pow(qq1*nu1,qq2*nu2,qq3*nu3,fit_ar[10],fit_ar[38],fit_ar[39],ell);}

                //if ((int)(fit_ar[12])==1){fog_pref = fog_B(qq1*nu1,qq2*nu2,qq3*nu3,sb);}
//                {fogk1 = fog_P_pow(qq1*nu1,sb,fit_ar[38],fit_ar[39],ell);
//                 fogk2 = fog_P_pow(qq2*nu2,sb,fit_ar[38],fit_ar[39],ell);
//                 fogk3 = fog_P_pow(qq3*nu3,sb,fit_ar[38],fit_ar[39],ell);}
//                {fog_pref = fog_B_atan(qq1*nu1,qq2*nu2,qq3*nu3,sb,fit_ar[38],fit_ar[39],ell);}
//                {fog_pref = fog_B_pow(qq1*nu1,qq2*nu2,qq3*nu3,sb,fit_ar[38],fit_ar[39],ell);}
            }

            z1_1 = z1spt(nu1,par_ar);
            z1_2 = z1spt(nu2,par_ar);
            z1_3 = z1spt(nu3,par_ar);

            z2_12 = fitz2spt(qq1,qq2,nu1,nu2, cos(tq12), par_ar, fit_ar,ell);
            z2_13 = fitz2spt(qq1,qq3,nu1,nu3, cos(tq13), par_ar, fit_ar,ell);
            z2_23 = fitz2spt(qq2,qq3,nu2,nu3, cos(tq23), par_ar, fit_ar,ell);

            temp12 = (2.0 * pq1 * pq2 * z1_1 * z1_2 * z2_12);
            temp13 = (2.0 * pq1 * pq3 * z1_1 * z1_3 * z2_13);
            temp23 = (2.0 * pq2 * pq3 * z1_2 * z1_3 * z2_23);
            
            mup_pr = new_mufacs(idx_mup,mu1,mu2,mu3);

            fval[0] = fog_pref * mup_pr * (temp12 + temp13 + temp23) / (pow(alpa*alpe*alpe,2.));
        }
        else {fval[0]=0.;}
        //if(idx_mup==1){printf("ka %.4f kb %.4f kc %.4f fog %.4f z213 %.4f z223 %.4f fval %.4f\n",qq1,qq2,qq3,fog_pref,z2_13,z2_23,fval[0]/(4. * M_PI));}
    }
    return 0;
}

extern void ext_bk_mqh(double **kh_ar, double *log_kh, double *log_pk, double *par_ar, double *fit_ar,
                       int kh_dim, int tr_num, int idx_mup, double **bk_mqh_ar)
{
    /// interpolator of the pk initialisation ///
    gsl_interp_accel *acc = gsl_interp_accel_alloc ();
    gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, kh_dim);
    gsl_spline_init (spline, log_kh, log_pk, kh_dim);
    //printf("f1 %.4f f2 %.4f f3 %.4f f4 %.4f f5 %.4f f6 %.4f f7 %.4f g1 %.4f g2 %.4f g3 %.4f g4 %.4f g5 %.4f g6 %.4f g7 %.4f sb %.4f\n",fit_ar[0],fit_ar[1],fit_ar[2],fit_ar[3],fit_ar[4],fit_ar[5],fit_ar[6],fit_ar[7],fit_ar[8],fit_ar[9],fit_ar[10],fit_ar[11],fit_ar[12],fit_ar[13],par_ar[6]);
    double xmin[2] = {-1.0,0.0}, xmax[2] = {1.0,2.0*M_PI},triangle[3]={0.,0.,0.},val, err;
    double xin[2] = {-0.3,1.1};
    for (int j = 0; j < idx_mup; j ++)
    {
        for (int i = 0 ; i < tr_num; i ++)
        {
//            multipole on the integrating mu
//            if (j==0)          {triangle[0] = kh_ar[i][0], triangle[1] = kh_ar[i][1], triangle[2] = kh_ar[i][2];}
//            if (j==1 || j == 4){triangle[0] = kh_ar[i][0], triangle[1] = kh_ar[i][1], triangle[2] = kh_ar[i][2];}
//            if (j==2 || j == 5){triangle[0] = kh_ar[i][1], triangle[1] = kh_ar[i][2], triangle[2] = kh_ar[i][0];}
//            if (j==3 || j == 6){triangle[0] = kh_ar[i][2], triangle[1] = kh_ar[i][0], triangle[2] = kh_ar[i][1];}

            //multipole on the third side, not directly related to the integrated variables
            if (j==0 || j == 7)          {triangle[0] = kh_ar[i][0], triangle[1] = kh_ar[i][1], triangle[2] = kh_ar[i][2];}
            if (j==1 || j == 4){triangle[0] = kh_ar[i][1], triangle[1] = kh_ar[i][2], triangle[2] = kh_ar[i][0];}
            if (j==2 || j == 5){triangle[0] = kh_ar[i][2], triangle[1] = kh_ar[i][0], triangle[2] = kh_ar[i][1];}
            if (j==3 || j == 6){triangle[0] = kh_ar[i][0], triangle[1] = kh_ar[i][1], triangle[2] = kh_ar[i][2];}
      
            //if(j==0){printf("ka %.4f kb %.4f kc %.4f \n",triangle[0],triangle[1],triangle[2]);}
            if(j==1){double z2_12 = fitz2spt(triangle[2],triangle[0],0.5,1.1, cos(tab_f(triangle[2], triangle[0], triangle[1])), par_ar, fit_ar,1);}
//            struct str_bkcbrt params = {kh_ar[i], par_ar, j, acc, spline,acc_n,spline_n,knl,par_ar[4]};
            struct str_bkcbrt params = {triangle, par_ar,fit_ar, j, acc, spline};
            bk_cbrt_gal(1,xin,&params,2,&val);
            //printf("ka %.4f kb %.4f kc %.4f fog %.4 \n",triangle[0],triangle[1],triangle[2],fog_B(ka,kb,kc,par_ar[6]);
            hcubature(1,bk_cbrt_gal, &params, 2, xmin, xmax, 0, 0, 5e-3,ERROR_INDIVIDUAL, &val, &err);
            bk_mqh_ar[j][i] = (1./(4. * M_PI)) * val;
//            printf("computed element %d multipole %d\n",i,j);
        }
    }
    // free memory of the interpolator
    gsl_spline_free (spline);
    gsl_interp_accel_free (acc);
}



