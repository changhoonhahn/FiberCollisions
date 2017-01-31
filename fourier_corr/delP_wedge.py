'''


Calculate the effective window method's correction on P(k) wedges. 


'''
import pickle
import numpy as np
from scipy.special import legendre
from scipy.integrate import quad as quad_int
import fourier_corr as four_corr 


def Tnl(en, EN, ell): 
    ''' Calculate the coefficient used for estiamte powerspectrum wedges from multipoles
   
   parameters
   ----------
        en : (int)
            en-th wedge out of EN equally spaced wedges for mu in [0, 1]
        EN : (int)
            Specifies the number of equally spaced wedges
        ell : multipole
    '''
    # calculate mu1 and mu2 first
    mu_wid = 1./np.float(EN)   
    mu1 = mu_wid*np.float(en - 1)
    mu2 = mu1 + mu_wid

    tmp_Leg = legendre(ell) # set up legendre polynomial
    return 1./(mu2 - mu1) * quad_int(tmp_Leg, mu1, mu2)[0]


def delPk_wedge(en, EN): 
    '''Calculate del P(k) wedge given P(k) multipole dictionary.  
    '''
    # hardcoded effective window method parameters
    fs = 0.6 
    rc = 0.43
    fold = 10
    rebin = 20

    for ell in [0, 2]:
        coeff = Tnl(en, EN, ell)

        # Del P Correlated
        k_corr, delP_corr = four_corr.DelPcorr_pkmu(ell, fs=fs, rc=rc, 
                fold=fold, rebin=rebin, dqdmu=True, dmudq=False)
       
        # Del P Uncorrelated
        k_uncorr, delP_uncorr = four_corr.DelPuncorr(ell, fs=fs, rc=rc, k_arr=k_corr)

        if ell == 0: 
            wedge = coeff * (delP_corr + delP_uncorr) 
        else: 
            wedge += coeff * (delP_corr + delP_uncorr) 

    return [k_corr, wedge]


def save_delPk_wedge(): 
    ''' Save delPk values to file  
    '''
    dP_dict = {} 
    for i in range(1,4): 
        dP_dict['delP_wedge_'+str(i)+'_3'] = delPk_wedge(i, 3)
    
    pickle.dump(dP_dict, open('/mount/riachuelo1/hahn/power/Nseries/Box/delP_3wedge_ell02.p', 'wb'))
    return None 


if __name__=='__main__': 
    #print delPk_wedge(1, 3)
    save_delPk_wedge()
