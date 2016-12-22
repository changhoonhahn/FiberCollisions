'''

P(k, mu) for Nseries box

'''
import time
import pickle
import numpy as np
import matplotlib as mpl 
import matplotlib.pyplot as plt
from scipy.io import FortranFile
from scipy.special import legendre
from scipy.special import j1 as J1
from scipy.integrate import quad as quad_int 
from scipy.integrate import dblquad as dblquad_int 
from scipy.interpolate import LinearNDInterpolator as LNDI

from Spectrum.interruptible_pool import InterruptiblePool as Pewl

import codif

from ChangTools.plotting import prettyplot 
from ChangTools.plotting import prettycolors 

def delPcorr_pkmu_file(ell, fs=0.7, rc=0.43, Ngrid=3600, fold=5, rebin=20, 
        dqdmu=True, dmudq=False, ktrust=None):
    ''' Pickle file name for Del P^corr from P(k, mu)
    '''
    if fold: 
        fold_str = ''.join(['.', str(fold), 'fold.', str(rebin), 'rebin'])
    else: 
        fold_str = ''

    if ktrust is None: 
        ktrust_str = ''
    else: 
        ktrust_str = ''.join(['.ktrust', str(round(ktrust,2))])

    if dqdmu and not dmudq: 
        integ_str = '.dqdmu'
    elif dmudq and not dqdmu: 
        integ_str = '.dmudq'
    else: 
        raise ValueError

    pickle_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'delPcorr_pkmu', 
        integ_str, 
        '.l', str(ell), 
        ktrust_str, 
        '.fs', str(round(fs,2)), 
        '.rc', str(round(rc,2)), 
        '.Ngrid', str(Ngrid), 
        fold_str, 
        '.p'])
    return pickle_file 

def delPcorr_pkmu_save(ell, fs=0.7, rc=0.43, Ngrid=3600, Nthreads=1, fold=5, rebin=20, 
        dqdmu=True, dmudq=False, ktrust=None):
    ''' Pickle Dump calculated Del P^corr from P(k, mu)
    '''
    pickle_file = delPcorr_pkmu_file(ell, fs=fs, rc=rc, Ngrid=Ngrid, fold=fold, rebin=rebin, 
            dqdmu=dqdmu, dmudq=dmudq, ktrust=ktrust)
    print 'Saving ', pickle_file 

    result = delPcorr_pkmu(ell, fs=fs, rc=rc, Ngrid=Ngrid, Nthreads=Nthreads, 
            fold=fold, rebin=rebin, dqdmu=dqdmu, dmudq=dmudq, ktrust=ktrust)

    pickle.dump(result, open(pickle_file, 'wb'))
    
    if ktrust is None: 
        ktrust_str = ''
    else: 
        ktrust_str = ''.join(['.ktrust = ', str(round(ktrust,2))])
    
    notif_msg = ''.join([
        'delP', str(ell), 'corr P(k,mu) ', 
        ktrust_str, ' ', 
        'fs='+str(round(fs,2)), ' ', 
        'rc='+str(round(rc,2))+' finished'])
    codif.notif(subject=notif_msg)
    return None

def delPcorr_pkmu(ell, fs=0.7, rc=0.43, Ngrid=3600, Nthreads=1, fold=5, rebin=20, 
        dqdmu=True, dmudq=False, ktrust=None): 
    ''' Calculate Del P^corr from P(k, mu)

    Parameters
    ----------
    ell : int
        Multiple l value. l = 0 or 2. 

    fs : float  
        Survey fraction that suffers from fiber collisions. fs < 1.0. For CMASS
        tiling fs ~ 0.7 

    rc : float
        Fiber collision comoving scale. At median redshift of CMASS it's 0.43

    Ngrid : int
        Grid size. Ngrid = 360 or 3600. Use 360 for testing and 3600 for actual 
        calculations 

    Nthreads : int
        Multiprocessing is implemented for the calculation. Specify the number of 
        threads you want to use  
    
    fold : int
        Number of folds in the supplementary P(k, mu). 

    rebin : int 
        Number of k bins for k > 0.1. This is to reduce the noise in P(k, mu)

    ktrust : float (optional)
        If specified, Del P^corr will only be integrated from 0 to ktrust with respect 
        to q.
    '''
    if dqdmu and dmudq: 
        raise ValueError

    # this is without the (2pi)^3 normalization
    #alphap = -4.0 * np.pi**3 * fs * rc**2 * (2. * ell + 1.)/2. 
    alpha = -0.5 * fs * rc**2 # normalized
    
    k_arr, mu_arr, avgP, n_mock = averaged_Pkmu(Ngrid=Ngrid, fold=fold, rebin=rebin)
    mask = ~np.isnan(avgP)

    Pkmu_interp = LNDI(np.vstack([mu_arr[mask], k_arr[mask]]).T, avgP[mask]) 
    #Pkmu_interp = interp2d(k_arr[mask], mu_arr[mask], avgP[mask], kind='linear')

    delP = [] 
    N_mubin = 20
    N_kbin = len(k_arr)/N_mubin
    k_bin = k_arr.reshape((N_kbin, N_mubin))[:,0]
    mu_bin = mu_arr.reshape((N_kbin, N_mubin))[0]

    integ_kwargs = {
            'ell': ell, 
            'rc': rc, 
            'Pkmu_interp': Pkmu_interp, 
            'k_bin': k_bin, 
            'mu_bin': mu_bin
            }
    if ktrust is not None: 
        integ_kwargs['ktrust'] = ktrust
        k_range = k_bin[np.where(k_bin < ktrust)]
    else: 
        k_range = k_bin

    if Nthreads > 1:  
        pool = Pewl(processes=Nthreads)
        mapfn = pool.map
        arglist = [[kii, integ_kwargs] for kii in k_range]
        
        if dqdmu: 
            delP = mapfn(delPcorr_integrate_wrap_dqdmu, [arg for arg in arglist])
        elif dmudq:
            delP = mapfn(delPcorr_integrate_wrap_dmudq, [arg for arg in arglist])
        pool.close()
        pool.terminate()
    else: 
        for kii in k_range: 
            if dqdmu: 
                delP_i = delPcorr_integrate_wrap_dqdmu([kii, integ_kwargs])
            elif dmudq: 
                delP_i = delPcorr_integrate_wrap_dmudq([kii, integ_kwargs])
            delP.append(delP_i)

    delP = np.array(delP)
    delP *= alpha * (2.*ell + 1.)/2.
    return [k_range, delP]

def delPcorr_integrate_wrap_dqdmu(args): 
    ''' Wrapper for the Del P^{corr} integration. The integrals are computed 
    in the following order : int (int ... dq) dmu
    '''
    kii = args[0]
    kwargs = args[1]
    
    t_integ = time.time()
    delPkmu_integ = lambda qq, mumu: delPcorr_integrand_kmu_dqdmu(
            qq, 
            mumu, 
            kii, 
            l=kwargs['ell'], 
            rc=kwargs['rc'], 
            Pkmu_interp=kwargs['Pkmu_interp'], 
            k_bin=kwargs['k_bin'], 
            mu_bin=kwargs['mu_bin']
            )
    
    if 'ktrust' in kwargs.keys(): 
        q_upper = kwargs['ktrust']
    else: 
        q_upper = kwargs['k_bin'][-1] 

    delP_i = dblquad_int(
            delPkmu_integ, 
            -1.,                        # mu lower lim 
            1.,                         # mu upper lim
            lambda x: np.abs(kii * x),   # q lower lim
            lambda x: q_upper,        # q upper lim 
            epsrel=0.005
            )

    print 'Del P'+str(kwargs['ell'])+'^corr(k='+str(round(kii,3))+')=', delP_i[0], ' t = ', time.time() - t_integ
    return delP_i[0]

def delPcorr_integrate_wrap_dmudq(args): 
    ''' Wrapper for the Del P^{corr} integration. The integrals are computed 
    in the following order : int (int ... dmu) dq
    '''
    kii = args[0]
    kwargs = args[1]
    
    t_integ = time.time()
    delPkmu_integ = lambda qq, mumu: delPcorr_integrand_kmu_dmudq(
            qq, 
            mumu, 
            kii, 
            l=kwargs['ell'], 
            rc=kwargs['rc'], 
            Pkmu_interp=kwargs['Pkmu_interp'], 
            k_bin=kwargs['k_bin'], 
            mu_bin=kwargs['mu_bin']
            )
    if 'ktrust' in kwargs.keys(): 
        q_upper = kwargs['ktrust']
    else: 
        q_upper = kwargs['k_bin'][-1] 
    delP_i = dblquad_int(
            delPkmu_integ, 
            kwargs['k_bin'][0],          # q lower lim 
            q_upper,                    # q upper lim
            lambda x: np.max([-1.0, -x/kii]),     # mu lower lim
            lambda x: np.min([1.0, x/kii]),       # mu upper lim 
            epsrel=0.005
            )

    print 'Del P'+str(kwargs['ell'])+'^corr(k='+str(round(kii,3))+')=', delP_i[0], ' t = ', time.time() - t_integ
    return delP_i[0]

def delPcorr_integrand_kmu_dqdmu(q, mu, k, l=None, rc=None, Pkmu_interp=None, k_bin=None, 
        mu_bin=None): 
    ''' Del P^corr integrand calculated using P(k, mu) and dblquad: 

    Leg_l(mu) q P(q, (k mu)/q) W_1D(k rc sqrt(1-mu^2), q rc sqrt(1 - (k mu/q)^2))
    '''
    #integ_time = time.time()
    Leg_l = legendre(l)
    theta = k*mu/q  
    x = k * rc * np.sqrt(1. - mu**2)
    y = q * rc * np.sqrt(1. - theta**2)

    # if q is beyond the k bounds
    if q < k_bin[0]: 
        return 0.0
    elif q > k_bin[-1]:
        return 0.0
    # to prevent boundary issues
    if theta > mu_bin[-1] and theta <= 1.0: 
        theta = mu_bin[-1]
    elif theta < mu_bin[0] and theta >= -1.0: 
        theta = mu_bin[0]

    w1d = W_2d(y) + (J1(y)/(y**2) - J1(y)/(2*y)) * x**2

    integrand = Leg_l(mu) * q * Pkmu_interp(np.array([theta, q]))[0] * w1d 
    #print 'delPcorr integrand takes ', time.time()-integ_time
    return integrand

def delPcorr_integrand_kmu_dmudq(mu, q, k, l=None, rc=None, Pkmu_interp=None, k_bin=None, 
        mu_bin=None): 
    ''' Del P^corr integrand calculated using P(k, mu) and dblquad: 

    Leg_l(mu) q P(q, (k mu)/q) W_1D(k rc sqrt(1-mu^2), q rc sqrt(1 - (k mu/q)^2))
    '''
    #integ_time = time.time()
    Leg_l = legendre(l)
    theta = k*mu/q  
    x = k * rc * np.sqrt(1. - mu**2)
    y = q * rc * np.sqrt(1. - theta**2)

    # if q is beyond the k bounds
    if q < k_bin[0]: 
        return 0.0
    elif q > k_bin[-1]:
        return 0.0
    # to prevent boundary issues
    if theta > mu_bin[-1] and theta <= 1.0: 
        theta = mu_bin[-1]
    elif theta < mu_bin[0] and theta >= -1.0: 
        theta = mu_bin[0]

    w1d = W_2d(y) + (J1(y)/(y**2) - J1(y)/(2*y)) * x**2

    integrand = Leg_l(mu) * q * Pkmu_interp(np.array([theta, q]))[0] * w1d 
    #print 'delPcorr integrand takes ', time.time()-integ_time
    return integrand


def W_2d(x): 
    '''
    W_2d in FC correction equation 

    W_2d = (2 J_1(x)/ x)
    '''
    return 2.*J1(x)/x

def averaged_Pkmu(Ngrid=3600, fold=5, rebin=20):
    ''' Given Ngrid, return averaged and symmetric P(k, mu) from Nseries 
    Box. Ouput is in the form [k_arr, mu_arr, sym_avgP]. All the arrays 
    have dimensions N_kbin * N_mubin. 

    Parameters
    ----------
    Ngrid : int
        Int that specifies the number of grids. Only 360 and 3600
    fold : int
        Int that specifies the folded P(k, mu) supplemented at high k
    rebin : int
        Int that specifies the number of k bins at k > 0.1. This is because the 
        folded P(k, mu) is very noisy.  
    '''
    if Ngrid == 3600:
        n_mock = 7
        n_mock_list = range(1, n_mock+1)
    elif Ngrid == 360: 
        n_mock = 1
        n_mock_list = [7] 
    else: 
        raise NotImplementedError
    
    for imock in n_mock_list: 
        k_i, mu_i, co_i, avgP_i = nseriesbox_pkmu(imock, Ngrid=Ngrid, ordered=False, 
                fold=fold, rebin=rebin)
        try: 
            nmock_avgP_i = np.repeat(1, len(avgP_i))
            nmock_avgP_i[np.isnan(avgP_i)] = 0
            nmock_avgP += nmock_avgP_i
            
            k_arr += k_i 
            mu_arr += mu_i 
            co_sum += co_i
            avgP_sum += avgP_i
            coP_sum += co_i * avgP_i

        except UnboundLocalError: 
            nmock_avgP = np.repeat(1, len(avgP_i))
            nmock_avgP[np.isnan(avgP_i)] = 0

            k_arr = k_i 
            mu_arr = mu_i 
            co_sum = co_i 
            avgP_i[np.isnan(avgP_i)] = 0
            avgP_sum = avgP_i
            coP_sum = co_i * avgP_i
    
    k_arr /= np.float(n_mock)
    mu_arr /= np.float(n_mock)
    allnan = np.where(nmock_avgP == 0)
    avgP = avgP_sum / nmock_avgP.astype(float)
    avgP[allnan] = np.NaN
    
    N_mubin = 20 
    N_kbin = len(k_arr)/N_mubin

    sym_avgP = avgP
    for ik in xrange(N_kbin): 
        for imu in xrange(N_mubin/2):
            ii = N_mubin * ik + imu
            sym_ii = N_mubin * ik + (19 - imu)
            
            if co_sum[ii] + co_sum[sym_ii] == 0: 
                sym_avgP_ii = np.NaN
            elif co_sum[ii] == 0 and co_sum[sym_ii] != 0: 
                sym_avgP_ii = coP_sum[sym_ii]/co_sum[sym_ii]
            elif co_sum[ii] != 0 and co_sum[sym_ii] == 0: 
                sym_avgP_ii = coP_sum[ii]/co_sum[ii]
            else:
                sym_avgP_ii = (coP_sum[ii] + coP_sum[sym_ii])/(co_sum[ii] + co_sum[sym_ii])

            sym_avgP[ii] = sym_avgP_ii
            sym_avgP[sym_ii] = sym_avgP_ii
    
    return [k_arr, mu_arr, sym_avgP, n_mock]

def nseriesbox_pkmu(imock, Ngrid=3600, ordered=True, fold=5, rebin=20):
    ''' Read k, mu, counts, and P(k,mu) from Nseries Box mocks and then
    reshape the arrays into a logical matrix: [k, mu, P(k,mu)]. If ordered is False, 
    then the k, mu, counts, and P(k,mu) are NOT reshaped. 

    Parameters
    ----------
    imock : int 
        Int that specifies box number (1 - 7) 
    Ngrid : int
        Int that specifies Ngrid. Ngrid can only be 360, or 3600. 
        360 for testing purposes.
    ordered : bool
        If True, return ordered array k_bin, mu_bin, and P(k, mu) 2D array. 
        If False, retrin 1 dimensional array. 
    fold : int
        Number of folds
    rebin : int
        Number of k bins past k > 0.1 once folded P(k, mu) are included. This
        is to make the P(k, mu) values less noisy. 
    '''
    # read k, mu, counts, P(k,mu)
    fname = ''.join(['/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'pkmu', str(Ngrid), 'z_BoxN', str(imock), '.dat'])
    k, mu, co, avgP = np.loadtxt(fname, unpack=True, usecols=[0,1,-2,-1])
    avgP *= (2*np.pi)**3
    k_max = round(k.max(), 3)

    # read k, mj counts from folded P(k,mu) 
    if fold: 
        fold_file = ''.join(['/mount/riachuelo1/hahn/power/Nseries/Box/', 
            'pkmu', str(Ngrid), 'z_BoxN', str(imock), '.', str(fold), 'xfold.dat'])
        k_f, mu_f, co_f, avgP_f = np.loadtxt(fold_file, unpack=True, usecols=[0,1,-2,-1])
        avgP_f *= (2*np.pi*np.float(fold))**3
        
        beyond_kmax = np.where(k_f > k_max)

        k = np.concatenate([k, k_f[beyond_kmax]])
        mu = np.concatenate([mu, mu_f[beyond_kmax]])
        co = np.concatenate([co, co_f[beyond_kmax]])
        avgP = np.concatenate([avgP, avgP_f[beyond_kmax]])
    
        N_mubin = 20
        if rebin: 
            norm_bin = np.where(k < 10.**0.)    # normal binning
            k_norm = k[norm_bin]
            mu_norm = mu[norm_bin]
            co_norm = co[norm_bin]
            avgP_norm = avgP[norm_bin]

            k_logbin = np.logspace(0.0, 2.0, rebin)
            for i_k in range(1, len(k_logbin)-1):
                k_logbin_lo = k_logbin[i_k-1]
                k_logbin_hi = k_logbin[i_k+1]

                within = np.where((k >= k_logbin_lo) & (k < k_logbin_hi))
                if len(within[0]) == 0: 
                    continue
                #print k_logbin_lo, k_logbin_hi, len(within[0])
                
                bin_k = k[within]
                bin_P = avgP[within]
                bin_co = co[within]
                bin_mu = mu[within]

                bin_k_2d = bin_k.reshape((len(bin_k)/N_mubin, N_mubin))
                bin_P_2d = bin_P.reshape((len(bin_P)/N_mubin, N_mubin))
                bin_co_2d = bin_co.reshape((len(bin_co)/N_mubin, N_mubin))
                bin_mu_2d = bin_mu.reshape((len(bin_mu)/N_mubin, N_mubin))

                sum_co_2d = np.sum(bin_co_2d, axis=0)
                avg_k_2d = np.sum(bin_k_2d * bin_co_2d, axis=0) / sum_co_2d
                avg_P_2d = np.sum(bin_P_2d * bin_co_2d, axis=0) / sum_co_2d
                avg_mu_2d = np.sum(bin_mu_2d * bin_co_2d, axis=0) / sum_co_2d 

                sum_co_1d = sum_co_2d.flatten()    
                avg_k_1d = avg_k_2d.flatten()
                avg_P_1d = avg_P_2d.flatten()
                avg_mu_1d = avg_mu_2d.flatten()

                k_norm = np.concatenate([k_norm, avg_k_1d])
                co_norm = np.concatenate([co_norm, sum_co_1d])
                mu_norm = np.concatenate([mu_norm, avg_mu_1d])
                avgP_norm = np.concatenate([avgP_norm, avg_P_1d])

            
            k = k_norm 
            mu = mu_norm
            co = co_norm 
            avgP = avgP_norm

    if not ordered:
        return k, mu, co, avgP
    else: 
        N_mubin = 20
        N_kbin = len(k)/N_mubin
        
        k_bin = k.reshape((N_kbin, N_mubin))[:,0]
        mu_bin = mu.reshape((N_kbin, N_mubin))[0]
        co = co.reshape((N_kbin, N_mubin))
        avgP = avgP.reshape((N_kbin, N_mubin))

        return k_bin, mu_bin, co, avgP

# --- Plotting --- 

def plot_nseriesbox_p_k_mu(imock, plot='2D', Ngrid=3600):
    ''' Plot P(k,mu). Option for plotting the monopole is included.
    '''
    # read k, mu, counts, P(k,mu)
    k_bin, mu_bin, co, avgP = nseriesbox_pkmu(imock, Ngrid=Ngrid)
    
    prettyplot()
    pretty_color = prettycolors()
    fig = plt.figure(figsize=(14,10))
    sub = fig.add_subplot(111)
    if plot == '2D': 
        colormap = plt.cm.afmhot
        cont = sub.pcolormesh(np.arange(-1.0, 1.1, 0.1), k_bin, avgP, cmap=colormap)
        plt.colorbar(cont)
        # y-axis
        sub.set_ylabel(r'$\mathtt{k \;(\mathrm{Mpc}/h})$', fontsize=20)
        # x-axis
        sub.set_xlabel(r'$\mu$', fontsize=20)
    elif plot == 'pk': 
        proj_p0k = [] 
        proj_p2k = [] 
        Leg2 = 2.5 * (3.0 * mu_bin**2 - 1.0)    # L_2(mu)
        for i in range(N_kbin): 
            pkmu_k = avgP[i]
            co_k = co[i]
            
            tot_p0k = np.sum(co_k[~np.isnan(pkmu_k)] * pkmu_k[~np.isnan(pkmu_k)])
            tot_count = np.sum(co_k[~np.isnan(pkmu_k)])
            proj_p0k.append(tot_p0k/tot_count)
            
            tot_p2k = np.sum(co_k[~np.isnan(pkmu_k)] * pkmu_k[~np.isnan(pkmu_k)] * Leg2[~np.isnan(pkmu_k)])
            tot_count = np.sum(co_k[~np.isnan(pkmu_k)])
            proj_p2k.append(tot_p2k/tot_count)
            
        sub.plot(k_bin, proj_p0k, lw=4, c=pretty_color[0], label='Monopole')
        sub.plot(k_bin, proj_p2k, lw=4, c=pretty_color[2], label='Quadrupole')
        k, p0k, p2k = np.loadtxt('/home/users/hahn/power360z_BoxN7.dat', unpack=True, usecols=[0,-1,2])
        sub.plot(k, p0k * (2*np.pi)**3, ls='--', c='k')
        sub.plot(k, p2k * (2*np.pi)**3, ls='--', c='k')

        # x-axis
        sub.set_xlabel(r'$\mathtt{k \;(\mathrm{Mpc}/h})$', fontsize=25)
        sub.set_xlim([10**-3, 10**1])
        sub.set_xscale("log")
        #y-axis
        sub.set_ylabel(r"$\mathtt{P_{l}(k)}$", fontsize=25)
        sub.legend(loc='upper right')
    else:
        raise ValueError

    sub.set_yscale("log")
    fig.savefig(''.join(['figure/', 
        'pkmu_', str(Ngrid), 'z_BoxN', str(imock), '.', plot, '.png']), 
        bbox_inches='tight')
    plt.close()

def plot_interp_p_k_mu(Ngrid=3600, proj=False, fold=5, rebin=20):
    ''' Compare pcolormesh plots of P(k, mu) and interpolated P(k,mu). 
    P(k, mu) averaged from 7 Nseries box mocks 

    Parameters
    ----------
    Ngrid : int
        Ngrid = 3600 or 360. 360 only for testing purposes. 
    proj : bool
        If True, plot P(k, mu) for a number of mu values. Else plot the 
        contour plot. 
    fold : int
        Number of folds in the supplementary small scale P(k, mu). Extends
        k_max by a factor of fold. 
    rebin : int
        Number of k bins for k > 0.1 once folded P(k, mu) is included. This 
        is to reduce the amount of noise. 
    '''
    k_arr, mu_arr, sym_avgP, n_mock = averaged_Pkmu(Ngrid=Ngrid, fold=fold, rebin=rebin)
    
    N_mubin = 20 
    N_kbin = len(k_arr)/N_mubin

    #mask = ~np.isnan(avgP)
    sym_mask = ~np.isnan(sym_avgP)
    
    # ordered P(k,mu)
    k_bin = k_arr.reshape((N_kbin, N_mubin))[:,0]
    mu_bin = mu_arr.reshape((N_kbin, N_mubin))[0]
    #avgPkmu = avgP.reshape((N_kbin, N_mubin))
    sym_avgPkmu = sym_avgP.reshape((N_kbin, N_mubin))
    
    # interpolation call function
    interp_P = LNDI(np.vstack([mu_arr[sym_mask], k_arr[sym_mask]]).T, sym_avgP[sym_mask]) 
    
    prettyplot()
    pretty_color = prettycolors()

    if proj: 
        fig = plt.figure()
        sub = fig.add_subplot(111)
        for ii in xrange(10): 
            i_rand = np.random.randint(0, 19)
            sub.plot(k_bin, sym_avgPkmu[:,i_rand], c='k') 
            k_bin_mid = 0.5 * (k_bin[1:] + k_bin[:-1])
            xi = np.vstack([np.repeat(0.5*(mu_bin[i_rand] + mu_bin[i_rand+1]), len(k_bin_mid)), k_bin_mid]).T
            if ii == 0: 
                label_str = 'Interpolated'
            else: 
                label_str = None 
            sub.plot(k_bin_mid, interp_P(xi), 'r', label=label_str)
            sub.plot(k_bin, sym_avgPkmu[:,i_rand+1], c='k') 
        # y-axis
        sub.set_yscale("log")
        sub.set_ylabel(r'$\mathtt{k \;(\mathrm{Mpc}/h})$', fontsize=30)
        # x-axis
        sub.set_xscale("log")
        sub.set_xlim([5.*10**-1, np.max(k_bin)])
        sub.set_xlabel(r'$\mu$', fontsize=30)
        sub.legend(loc='upper right') 
        proj_str = '.proj'
    else: 
        fig = plt.figure(figsize=(20,10))
        for i_sub in [1, 2]: 
            sub = fig.add_subplot(1, 2, i_sub)
            colormap = plt.cm.afmhot
            if i_sub == 1: 
                print 'Min[P(k, mu)] = ', sym_avgP[sym_mask].min()
                print 'Max[P(k, mu)] = ', sym_avgP[sym_mask].max()
                #vmin=sym_avgP[sym_mask].min(), 
                #vmax=sym_avgP[sym_mask].max(), 
                norm = mpl.colors.SymLogNorm(1000., 
                        vmin=-20., 
                        vmax=127675.,
                        clip=True)
                cont = sub.pcolormesh(np.arange(-1.0, 1.1, 0.1), k_bin, sym_avgPkmu, 
                        norm=norm, cmap=colormap)
                sub.set_ylabel(r'$\mathtt{k \;(\mathrm{Mpc}/h})$', fontsize=30)
                plt.colorbar(cont)
            else: 
                print 'P(k = ', k_bin[10], ', mu = ', mu_bin[10], ') = ', sym_avgPkmu[10,10]
                print 'P(k = ', k_bin[10], ', mu = ', mu_bin[10], ') Interpolated = ', interp_P([mu_bin[10], k_bin[10]])

                Pkmu_interp = interp_P(np.vstack([mu_arr, k_arr]).T)
                norm = mpl.colors.SymLogNorm(1000., 
                        vmin=-20., 
                        vmax=127675.,
                        clip=True)
                Pkmu_interp = Pkmu_interp.reshape((N_kbin, N_mubin))
                cont = sub.pcolormesh(np.arange(-1.0, 1.1, 0.1), k_bin, Pkmu_interp, 
                        norm=norm, cmap=colormap)
                plt.colorbar(cont)
            # y-axis
            sub.set_yscale("log")
            sub.set_ylim([np.min(k_bin), np.max(k_bin)])
            # x-axis
            sub.set_xlim([-1., 1.])
            sub.set_xlabel(r'$\mu$', fontsize=30)
        proj_str = ''

    if fold: 
        fold_str = ''.join(['.', str(fold), 'fold'])
    else: 
        fold_str = ''
    
    fig.savefig(''.join(['figure/', 
        'interp_pkmu_', str(Ngrid), 'z_BoxN', 
        '.', str(n_mock), 'mock', 
        fold_str, proj_str, '.png']), 
        bbox_inches='tight')
    plt.close()

def plot_delPcorr_integrand(ell, k, fs=0.7, rc=0.43, Ngrid=3600, fold=5, rebin=20):
    '''
    '''
    if Ngrid == 3600: 
        n_mock = 7
        n_mock_list = range(1, n_mock+1)
    elif Ngrid == 360:
        n_mock = 1
        n_mock_list = [7]

    for imock in n_mock_list: 
        k_i, mu_i, co_i, avgP_i = nseriesbox_pkmu(imock, Ngrid=Ngrid, ordered=False, 
                fold=fold, rebin=rebin)
        try: 
            k_arr += k_i 
            mu_arr += mu_i 
            avgP += avgP_i 
        except UnboundLocalError: 
            k_arr = k_i 
            mu_arr = mu_i 
            avgP = avgP_i 

    k_arr /= np.float(n_mock) 
    mu_arr /= np.float(n_mock) 
    avgP /= np.float(n_mock) 
    mask = ~np.isnan(avgP)
    Pkmu_interp = LNDI(np.vstack([mu_arr[mask], k_arr[mask]]).T, avgP[mask]) 
    
    N_mubin = 20
    N_kbin = len(k_arr)/N_mubin
    k_bin = k_arr.reshape((N_kbin, N_mubin))[:,0]
    mu_bin = mu_arr.reshape((N_kbin, N_mubin))[0]
    
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(10,7))
    sub = fig.add_subplot(111)

    q_arr = np.logspace(-3, 1, num=100)
    mu_arr = np.linspace(-0.95, 0.95, num=5)
    
    integ_q = []
    for i_mu, mu_i in enumerate(mu_arr): 
        integrand_q, q_vals = [], [] 
        for qi in q_arr: 
            if qi > np.abs(k * mu_i): 
                integrand_q.append(
                        delPcorr_integrand_kmu(
                            qi, mu_i, k, l=ell, rc=rc, Pkmu_interp=Pkmu_interp, k_bin=k_bin, 
                            mu_bin=mu_bin)
                        )
                q_vals.append(qi)
            else: 
                continue 
        integ_q.append([q_vals, integrand_q])
    for iint, int_q in enumerate(integ_q): 
        sub.plot(int_q[0], int_q[1], 
                lw=2, c=pretty_colors[iint], label=r'$\mu = '+str(round(mu_arr[iint],2))+'$')
    # x-axis
    sub.set_xscale('log')
    sub.set_xlim([10**-3, 10**2])
    sub.set_xlabel(r'$\mathtt{q}$', fontsize=25)
    # y-axis
    if ell == 0: 
        sub.set_ylim([0., 5000.])
    elif ell == 2: 
        sub.set_ylim([-2000., 2000.])

    sub.legend(loc='upper right') 
    
    fig_file = ''.join(['figure/', 
        'delP', str(ell), 'corr_pkmu_integrand', 
        '.k', str(round(k, 2)), 
        '.fs', str(round(fs, 2)),
        '.rc', str(round(rc, 2)), 
        '.Ngrid', str(Ngrid), 
        '.png'])

    fig.savefig(fig_file, bbox_inches='tight', dpi=200)
    plt.close()
    return None

def plot_delPcorr_mu_integrated(ell, fs=0.6, rc=0.43, Ngrid=3600, fold=5, rebin=20): 
    ''' In order to determine the effect of q_max = 4.3 cut-off or 
    potential underestimation of P(k, mu) for low k's plot the mu integrated
    integrand of Del P^corr for a range of k values. 

    If the effect of q_max cut-off is relatively constant for k < 0.1
    and more severe at low k's compared to high k's. The discrepancy in 
    Del P may be due to the cut-off. 
    '''
    if Ngrid == 3600: 
        n_mock = 7
        n_mock_list = range(1, n_mock+1)
    elif Ngrid == 360:
        n_mock = 1
        n_mock_list = [7]

    for imock in n_mock_list: 
        k_i, mu_i, co_i, avgP_i = nseriesbox_pkmu(imock, Ngrid=Ngrid, ordered=False, 
                fold=fold, rebin=rebin)
        try: 
            k_arr += k_i 
            mu_arr += mu_i 
            avgP += avgP_i 
        except UnboundLocalError: 
            k_arr = k_i 
            mu_arr = mu_i 
            avgP = avgP_i 

    k_arr /= np.float(n_mock) 
    mu_arr /= np.float(n_mock) 
    avgP /= np.float(n_mock) 
    mask = ~np.isnan(avgP)
    Pkmu_interp = LNDI(np.vstack([mu_arr[mask], k_arr[mask]]).T, avgP[mask]) 
    
    N_mubin = 20
    N_kbin = len(k_arr)/N_mubin
    k_bin = k_arr.reshape((N_kbin, N_mubin))[:,0]
    mu_bin = mu_arr.reshape((N_kbin, N_mubin))[0]
    
    prettyplot()
    pretty_colors = prettycolors()
    fig = plt.figure(figsize=(10,7))
    sub = fig.add_subplot(111)

    q_arr = np.logspace(-3, 1, num=100)     
    dq_integrand = np.zeros(len(q_arr)) 
    for i_k, k in enumerate([0.01, 0.05, 0.1, 0.3, 0.5, 0.8]): 
        for i_q, qi in enumerate(q_arr):  
            mu_integ = lambda mumu: delPcorr_integrand_kmu(
                    qi,                                 # q
                    mumu,                               # mu (the variable)
                    k,                                  # k 
                    l=ell,                              # ell 
                    rc=rc,                              # r 
                    Pkmu_interp=Pkmu_interp,            # interpolation function 
                    k_bin=k_bin, 
                    mu_bin=mu_bin
                    )
            lower_bound = np.max([-1., -1.*qi/k])
            upper_bound = np.min([1., qi/k])
            dq_integ = quad_int(mu_integ, lower_bound, upper_bound)
            dq_integrand[i_q] = dq_integ[0]

        sub.plot(q_arr, dq_integrand * 0.5 * (2.0 * np.float(ell) + 1.), 
                lw=2, c=pretty_colors[i_k], 
                label=r"$\mathtt{k = "+str(round(k,2))+"}$"
                )
    # x-axis
    sub.set_xscale('log')
    sub.set_xlim([10**-3, 10**2])
    sub.set_xlabel(r'$\mathtt{q}$', fontsize=25)
    # y-axis
    sub.set_ylabel(r"$\frac{2l + 1}{2} \int \mathcal{L}_l(\mu) q P(q, \theta) W_{1D} d\mu$", 
            fontsize=20)
    if ell == 0: 
        sub.vlines(k_bin[-1], 0., 4000., linestyle='--', linewidth=3, label=r"$\mathtt{q_{max}}$")
        sub.set_ylim([0.0, 2000.0])
    elif ell == 2: 
        sub.vlines(k_bin[-1], -800., 400., linestyle='--', linewidth=3, label=r"$\mathtt{q_{max}}$")
        sub.set_ylim([-2000.0, 1000.0])
    sub.legend(loc='upper right')

    if fold: 
        fold_str = ''.join(['.', str(fold), 'fold.', str(rebin), 'rebin'])
    else: 
        fold_str = ''
    
    fig_file = ''.join(['figure/', 
        'delP', str(ell), 'corr_pkmu', 
        '.1D_mu_integrated', 
        '.fs', str(round(fs, 2)),
        '.rc', str(round(rc, 2)), 
        '.Ngrid', str(Ngrid), 
        fold_str, 
        '.png'])
    fig.savefig(fig_file, bbox_inches='tight', dpi=150)
    plt.close()
    return None



"""
    DEPRECIATED BECAUSE THE 2D INTEGRAL IS COMPUTED DIRECTLY 
    def delPcorr_integrand(mu, k, rc=0.43, l=0, Pkmu_interp=None, k_bin=None, mu_bin=None):
        ''' Integrand of Del P^corr using P(k,mu) instead of the multipoles 

        Del P integrand = Leg_l(mu) * \int_kmu^inf q P(q, k*mu/q) W_2d(rc sqrt(q^2 - k^2mu^2))
        '''
        #int_time = time.time()
        Leg_l = legendre(l)
        qPqW2d_int = lambda qq: qPqW2d(
                qq,                                 # q 
                k,                                  # k  
                mu,                                 # mu 
                rc=rc,                              # rc
                Pkmu_interp=Pkmu_interp,        # P(k,mu) interpolation function
                k_bin=k_bin,                        # k bins
                mu_bin=mu_bin                       # mu bins
                )
        delPinteg = Leg_l(mu) * quad_int(qPqW2d_int, np.abs(k * mu), k_bin[-1], limit=10)
        #print 'DelP corr integrand ', time.time()-int_time
        return delPinteg[0]

    def qPqW2d(q, k, mu, rc=0.43, Pkmu_interp=None, k_bin=None, mu_bin=None): 
        '''
        q * P(q, theta) * W2D(q*rc*sqrt(1-theta^2))
        where theta = k*mu/q
        '''
        #qPqW2d_time = time.time()
        theta = k*mu/q  
        # if q is beyond the k bounds
        if q < k_bin[0]: 
            return 0.0
        elif q > k_bin[-1]:
            return 0.0
        # to prevent boundary issues
        if theta > mu_bin[-1] and theta <= 1.0: 
            theta = mu_bin[-1]
        elif theta < mu_bin[0] and theta >= -1.0: 
            theta = mu_bin[0]
        integrand = q * Pkmu_interp(np.array([theta,q])) * W_2d(q * rc * np.sqrt(1. - theta**2))
        #print 'qPqW2d ', time.time() - qPqW2d_time
        return integrand

    def read_fortran_fft(fname, N_bin): 
        ''' Read fortran FFT output 
        '''
        f = FortranFile(fname, 'r')
        record = f.read_record([
            ('Lx', '<i4'), 
            ('Ly', '<i4'), 
            ('Lz', '<i4'), 
            ('Npar', '<i4'), 
            ('akx', '<f4'), 
            ('aky', '<f4'), 
            ('akz', '<f4'), 
            ('phys_nyq', '<f4')
            ])
        Lx = record['Lx'][0]
        Ly = record['Ly'][0]
        Lz = record['Lz'][0]
        Npar = record['Npar'][0]
        akx = record['akx'][0]
        aky = record['aky'][0]
        akz = record['akz'][0]
        phys_nyq = record['phys_nyq'][0]

        print Lx, Ly, Lz, Npar, akx, aky, akz, phys_nyq
        dclr=f.read_reals(dtype=np.complex64)
        dclr=np.reshape(dclr,(Lx//2+1,Ly,Lz),order='F')

        mu_bins = np.arange(-1.0, 1.1, 0.1)
        mu_arr = np.arange(-0.95, 1.05, 0.1)
        N_mubin = len(mu_arr)

        co = np.zeros((N_bin, N_mubin))
        cok = np.zeros(N_bin)
        comu = np.zeros(N_mubin)
        
        avgk = np.zeros(N_bin)
        avgmu = np.zeros(N_mubin)
        avgP = np.zeros((N_bin, N_mubin))
        
        # k_z
        rkz = akz * (((np.arange(1, Lz+1) + Lz/2 - 2) % Lz) - Lz/2 + 1).astype(float)
        rkz_3d = np.repeat( np.repeat(rkz[np.newaxis, :], Ly, axis=0)[np.newaxis, :, :], Lx, axis=0)
        # k_y
        rky = aky * (((np.arange(1, Ly+1) + Ly/2 - 2) % Ly) - Ly/2 + 1).astype(float)
        rky_3d = np.repeat( np.repeat(rky[:, np.newaxis], Lz, axis=1)[np.newaxis, :, :], Lx, axis=0)
        # k_x 
        rkx = akx * (((np.arange(1, Lx+1) + Lx/2 - 2) % Lx) - Lx/2 + 1).astype(float)
        rkx_3d = np.repeat( np.repeat(rkx[:, np.newaxis], Lz, axis=1)[:, np.newaxis, :], Ly, axis=1)
        # k 
        rk_3d = np.sqrt(rkx_3d**2 + rky_3d**2 + rkz_3d**2)
        rk_notzero = np.where(rk_3d > 0.)
        imk_3d = (np.rint(np.float(N_bin) * rk_3d / phys_nyq)).astype(int)
        # mu
        mu_3d = rkz_3d / rk_3d
        imu_3d = np.rint((mu_3d + 1.05)/0.1).astype(int)

        ct = np.zeros((Lx, Ly, Lz), dtype=np.complex64)
        ct[:Lx//2+1, :, :] = dclr
        ct[:Lx//2:-1,Ngrid:0:-1,Ngrid:0:-1]=np.conj(phi[1:Ngrid//2,1:Ngrid,1:Ngrid])



        count = np.zeros((N_bin, N_mubin))
        avgPk = np.zeros((N_bin, N_mubin))
        for imk in range(1, N_bin+1): 
            for imu in range(1, N_mubin=1): 
                k_mu_bin = np.where((imk_3d == imk) & (imu_3d == imu))
                count = len(k_mu_bin[0])
                avgPk
        return None
"""


if __name__=='__main__': 
    #plot_nseriesbox_p_k_mu(1, plot='2D)
    #plot_interp_p_k_mu(Ngrid=3600, proj=False, fold=5)
    #plot_interp_p_k_mu(Ngrid=3600, proj=True, fold=5)
    #for ell in [0, 2]: 
    #    plot_delPcorr_mu_integrated(ell, fs=0.6, rc=0.43, Ngrid=3600, fold=10, rebin=20)
    #    plot_delPcorr_mu_integrated(ell, fs=0.6, rc=0.43, Ngrid=3600, fold=False, rebin=False)

    for ell in [2]: 
        delPcorr_pkmu_save(ell, fs=0.6, rc=0.43, Ngrid=3600, fold=10, rebin=20, 
                dmudq=False, dqdmu=True, ktrust=None, Nthreads=20)
