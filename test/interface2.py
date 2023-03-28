import numpy as np
import ctypes
from ctypes import *
from numpy.ctypeslib import ndpointer
from scipy.interpolate import interp1d
from scipy.signal import savgol_filter

"define a pointer for 1D arrays"
_doublep  = ndpointer(ctypes.c_double, flags='C_CONTIGUOUS')
"define a pointer for 1D arrays INT "
_intp  = ndpointer(ctypes.c_int, flags='C_CONTIGUOUS')
"define a pointer for 2D arrays"
_doublepp = ndpointer(dtype=np.uintp, ndim=1, flags='C')
"function to convert 2D array into a proper format for C"
def c_2d_inp(x):
    return (x.__array_interface__['data'][0]
            + np.arange(x.shape[0]) * x.strides[0]).astype(np.uintp)
path_lib = '/home/sergi/Geo-FPT/C_lib2/'
clib     = CDLL(path_lib + 'clib2.so')

""" BK MULTIPOLES"""
clib.ext_bk_mp.restype  = None
clib.ext_bk_mp.argtypes = [_doublepp, _doublepp,_doublepp,_doublepp, 
                           _doublep, _doublep, _doublep, _doublep,_doublep,_doublep,_doublep, _doublep, _doublep, c_int,  c_int, c_int, 
                           c_int,c_int,c_int,_doublep, c_double,c_double,c_int,c_double,c_double,c_int,c_int]


def bk_multip(tr, tr2,tr3,tr4,log_kli, log_pkli, log_kp, log_pk, cosm_par, af,ag,af_mup,ag_mup,D_z,amx,mod,sig_fog,sig_fog2,multip,mat):
    bk_out = np.zeros(tr.shape[0]+tr2.shape[0]+tr3.shape[0]+tr4.shape[0], dtype='float')
    kden = log_kp
    pk_den = log_pk
    knl, log_kh_n, n_smooth = eff_ingr(log_kli,log_pkli)
    clib.ext_bk_mp(c_2d_inp(tr), c_2d_inp(tr2),c_2d_inp(tr3),c_2d_inp(tr4),kden, pk_den, 
                   log_kh_n, n_smooth, cosm_par, af, ag,af_mup,ag_mup,kden.size, log_kh_n.size, tr.shape[0],tr2.shape[0],
                   tr3.shape[0],tr4.shape[0], bk_out, knl,D_z, mod,sig_fog,sig_fog2,multip,mat)
        
    return  bk_out

def multips(tr,tr2,tr3,tr4,k_lin,pk_lin, k_p, pk, cosm_par,x, x2, D_z,redshifts,full,multip,mod,mat):
    amx = 0.01
    if full ==1:
        if mod == 0: af,ag = x,x[3:-1]
        if mod == 1 or mod==9: af,ag = x,x #x[9:-1]
        if mod == 3: af,ag = x[:5],np.zeros(5) #,x[5:-1]
        if mod == 2: af,ag = x,np.zeros(8)
        if mod == 5: af,ag = x,np.zeros(8)
        if mod == 10: af,ag= x,np.zeros(6)
        af_mup,ag_mup=af,ag

    if full ==0:
        ag = x
        if mod==3:
#             if redshifts==1: af = np.array([ 1.00321199,  0.00778061,  0.04035794, -0.07043621,  0.01390807, 4.71829494])
#             if redshifts==3: af = np.array([ 1.15615701,  0.01335707, -0.00582238, -0.08951067,  0.01584383, 3.98633874])
#             if redshifts==4: af = np.array([ 1.29787528,  0.01154473, -0.06984416, -0.09486691,  0.01642198, 2.61686876])
            if redshifts==1: 
                af = np.array([ 1.01929699, -0.00595009,  0.02563233, -0.07043232,  0.01530769,5.06958099])
                af = np.array([ 1.01929699, -0.00595009,  0.02563233, -0.07043232,  0.01530769,5.06958099])
            if redshifts==3: af = np.array([ 1.03017161e+00, -4.73338583e-05,  2.36928681e-02, -6.90831429e-02,1.25805008e-02, 4.23352743e+00])
            if redshifts==4: af = np.array([ 1.06657406e+00, -7.49512475e-03,  5.18226238e-04, -6.29847783e-02,1.05734684e-02,  2.98461773e+00])
    af_mup,ag_mup=af,ag
    sig_fog = 0.
    if mod==3:
        if multip>=0: sig_fog = x[-1]
#        if full==1: sig_fog = cosm_par[4]
        if full==0: sig_fog = af[-1]
        sig_fog = cosm_par[4]
    sig_fog2 =sig_fog
    bk_multip = bk_multip(tr, tr2,tr3,tr4,np.log10(k_lin), np.log10(pk_lin), np.log10(k_p), np.log10(pk), cosm_par, af, ag, af,ag,D_z, amx,mod,sig_fog,sig_fog2,multip,mat)

    ind,ind2,ind3 = tr.shape[0],tr2.shape[0],tr3.shape[0]

    bk0   = bk_multip[:ind]
    bk200 = bk_multip[ind:(ind+ind2)]
    bk020 = bk_multip[(ind+ind2):(ind+ind2+ind3)]
    bk002 = bk_multip[(ind+ind2+ind3):]

    return bk0, bk200, bk020, bk002




#clib.ext_Preal.restype  = None
#clib.ext_Preal.argtypes = [_doublepp,_doublep,_doublep,c_double,c_int]

#def Pdd(theory,kin,cosm_par):
#    Pout = np.zeros(kin.size)
#    sigma8_sc = (cosm_par[0]/theory[-1,-1])**2
#    Num = kin.size
    #print(kin,sigma8_sc,Num)




clib.ext_Preal.restype  = None
clib.ext_Preal.argtypes = [_doublepp,_doublep,_doublep,c_double,c_double,c_double,c_int]

def Pdd(theory,kin,cosm_par):
    Pout = np.zeros(kin.size)
    alpa,alpe = 1.,1.#cosm_par[5],cosm_par[6]
    sigma8_sc = (cosm_par[0]/theory[-1,-1])**2
    Num = kin.size
    clib.ext_Preal(c_2d_inp(theory),kin,Pout,alpa,alpe,sigma8_sc,Num)
    return Pout
def Pdd_dav(theory,kin,cosm_par):
    Pout = np.zeros(kin.size)
    alpa,alpe = 1.,1.#cosm_par[4],cosm_par[5]
    sigma8_sc = (cosm_par[2]/theory[-1,-1])**2
    Num = kin.size
    clib.ext_Preal(c_2d_inp(theory),kin,Pout,alpa,alpe,sigma8_sc,Num)
    return Pout

clib.ext_Prsd.restype  = None
clib.ext_Prsd.argtypes = [_doublepp,_doublep,c_int,_doublep,c_double,c_double,c_double,c_int,c_double,c_double,c_double ,c_double,c_double,c_double,c_double,c_int,c_double]

def P024(theory,kin,cosm_par,mup):
    bs2 = -4./7.*(cosm_par[1]-1.)
    b3nl= 32./315.*(cosm_par[1]-1.)
    A = 0.#cosm_par[-2]
    Pnoise = 0. # 2393.84 
    Pout = np.zeros(kin.size)
    #print(Pout.shape)
    b1,b2,f,sig_P = cosm_par[1],cosm_par[3],cosm_par[2],cosm_par[-1]
    alpa,alpe = cosm_par[5],cosm_par[6]
    sigma8_sc = (cosm_par[0]/theory[-1,-1])**2
    Num = Pout.size
    kin_dim = kin.size
    
    clib.ext_Prsd(c_2d_inp(theory),kin,kin_dim,Pout,alpa,alpe,sigma8_sc,Num,
                  b1,b2,bs2,b3nl,A,Pnoise,f,mup,sig_P)
    return Pout






def sqrt_area(l):
    pos = 2.0 * ((l[0]*l[1])**2 +(l[0]*l[2])**2 + (l[1]*l[2])**2)
    neg = l[0]**4 + l[1]**4 + l[2]**4
    return np.sqrt(pos - neg)/4.

def eff_ingr(log_kh, log_pk):
    "compute n_smooth and knl to pass to the effective bispectrum formula"
    log_kh_n, n_smooth = smooth_n(log_kh, log_pk, 0.005)
    knl = 10 ** log_kh[(np.abs(((10 ** log_pk) * (10 ** log_kh) ** 3) / (2.0 * np.pi ** 2) - 1.0)).argmin()]
    return knl, log_kh_n, n_smooth


"""""""""""""""""""""""""""""""""""""""
Return smoothed slope of the the linear PK
"""""""""""""""""""""""""""""""""""""""
" four point interpolation"
def four_p_der(x, y, h):

    f = interp1d(x,y)

    # print(x.max(),x[-10]+2.*h)

    der = f(x[10:-10] -2.0 * h) + 8.0 * f(x[10:-10] + h) - (8.0 * f(x[10:-10] - h) + f(x[10:-10] +2.0 * h))

    return der/(12. * h)


def smooth_n(log_kh, log_pk, h):
    n0 = four_p_der(log_kh, log_pk, h)
    x_n0 = log_kh[10:-10]

    n1 = four_p_der(x_n0, n0, h)
    x_n1 = x_n0[10:-10]

    n2 = four_p_der(x_n1, n1, h)
    x_n2 = x_n1[10:-10]

    "find zeros second derivative"
    asign = np.sign(n2)
    sz = asign == 0
    while sz.any():
        asign[sz] = np.roll(asign, 1)[sz]
        sz = asign == 0

    signchange = ((np.roll(asign, 1) - asign) != 0).astype(bool)

    signchange[0] = 0

    x_n2_red = x_n2[signchange]

    ind_a = np.argwhere(x_n0 == x_n2_red[0])[0][0]
    ind_c = np.argwhere(x_n0 == x_n2_red[-1])[0][0]

    b_part_x = x_n0[ind_a:ind_c + 1]
    b_part_y = n0[ind_a:ind_c + 1]

    "find all local minima and maxima"
    indexes_min = np.r_[True, b_part_y[1:] < b_part_y[:-1]] & np.r_[b_part_y[:-1] < b_part_y[1:], True]
    indexes_max = np.r_[True, b_part_y[1:] > b_part_y[:-1]] & np.r_[b_part_y[:-1] > b_part_y[1:], True]

    a_bmin_c_x = np.r_[x_n0[:ind_a], b_part_x[indexes_min], x_n0[ind_c + 1:]]
    a_bmin_c_y = np.r_[n0[:ind_a], b_part_y[indexes_min], n0[ind_c + 1:]]

    f_min = interp1d(a_bmin_c_x, a_bmin_c_y)

    a_bmax_c_x = np.r_[x_n0[:ind_a], b_part_x[indexes_max], x_n0[ind_c + 1:]]
    a_bmax_c_y = np.r_[n0[:ind_a], b_part_y[indexes_max], n0[ind_c + 1:]]

    f_max = interp1d(a_bmax_c_x, a_bmax_c_y)

    n_min = f_min(x_n0)
    n_max = f_max(x_n0)

    n_smooth = 0.5 * (n_min + n_max)

    # print( int(np.rint(log_kh.size /3.0)))

    # n_smooth2 = savgol_filter(n_smooth, np.rint(log_kh.size /3.0), 3)
    n_smooth2 = savgol_filter(n_smooth, int(np.rint(log_kh.size /3.0)), 3)

    return x_n0, n_smooth2
