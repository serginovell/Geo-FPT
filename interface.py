import numpy as np
import ctypes
from ctypes import *
from numpy.ctypeslib import ndpointer

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
path_lib = './C_lib/' 
clib     = CDLL(path_lib + 'clib.so')

""" BK MULTIPOLES"""
clib.ext_bk_mp.restype  = None
clib.ext_bk_mp.argtypes = [_doublepp, _doublepp,_doublepp,_doublepp,_doublep, _doublep,_doublep, c_double,c_int,c_int,c_int,c_int,c_int,c_int, _doublep]


def bk_multip(tr, tr2,tr3,tr4, kp, pk, cosm_par, redshift,fit_full=1):
    bk_out = np.zeros(tr.shape[0]+tr2.shape[0]+tr3.shape[0]+tr4.shape[0], dtype='float')

    clib.ext_bk_mp(c_2d_inp(tr), c_2d_inp(tr2),c_2d_inp(tr3),c_2d_inp(tr4),np.log10(kp), np.log10(pk), 
                   cosm_par, redshift, fit_full, kp.size, tr.shape[0],tr2.shape[0],
                   tr3.shape[0],tr4.shape[0], bk_out)
    
    ind,ind2,ind3 = tr.shape[0],tr2.shape[0],tr3.shape[0]

    bk0   = bk_out[:ind]
    bk200 = bk_out[ind:(ind+ind2)]
    bk020 = bk_out[(ind+ind2):(ind+ind2+ind3)]
    bk002 = bk_out[(ind+ind2+ind3):]

        
    return  bk0, bk200, bk020, bk002





#clib.ext_Preal.restype  = None
#clib.ext_Preal.argtypes = [_doublepp,_doublep,_doublep,c_double,c_int]

#def Pdd(theory,kin,cosm_par):
#    Pout = np.zeros(kin.size)
#    sigma8_sc = (cosm_par[0]/theory[-1,-1])**2
#    Num = kin.size
    #print(kin,sigma8_sc,Num)
#    clib.ext_Preal(c_2d_inp(theory),kin,Pout,sigma8_sc,Num)
#    return Pout

clib.ext_Preal.restype  = None
clib.ext_Preal.argtypes = [_doublepp,_doublep,_doublep,c_double,c_double,c_double,c_int]

def Pdd(theory,kin,cosm_par):
    Pout = np.zeros(kin.size)
    alpa,alpe = 1.,1.
    sigma8_sc = (cosm_par[0]/theory[-1,-1])**2
    Num = kin.size
    clib.ext_Preal(c_2d_inp(theory),kin,Pout,alpa,alpe,sigma8_sc,Num)
    return Pout

clib.ext_Prsd.restype  = None
clib.ext_Prsd.argtypes = [_doublepp,_doublep,c_int,_doublep,c_double,c_double,c_double,c_int,c_double,c_double,c_double ,c_double,c_double,c_int,c_double]

def P024(theory,kin,cosm_par,mup):
    bs2 = -4./7.*(cosm_par[4]-1.)
    b3nl= 32./315.*(cosm_par[4]-1.)
    b1,b2,f,sig_P = cosm_par[4],cosm_par[5],cosm_par[1],cosm_par[7]
    alpa,alpe = cosm_par[2],cosm_par[3]
    sigma8_sc = (cosm_par[0]/theory[-1,-1])**2

    
    Pout = np.zeros(kin.size)
    Num = Pout.size
    kin_dim = kin.size
    
    clib.ext_Prsd(c_2d_inp(theory),kin,kin_dim,Pout,alpa,alpe,sigma8_sc,Num,
                  b1,b2,bs2,b3nl,f,mup,sig_P)

    return Pout

