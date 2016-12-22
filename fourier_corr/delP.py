'''

Module for TopHat convolution for Del P

'''
__author__=['ChangHoon Hahn']

import time 
import numpy as np
import matplotlib.pyplot as plt
import codif

import mpfit
import pickle 
import pk_extrap
import fourier_corr

from Spectrum.interruptible_pool import InterruptiblePool as Pewl

from scipy.special import j1 as J1
from scipy.interpolate import interp1d
from scipy.integrate import quad as quad_int 
from scipy.integrate import nquad 


def delPuncorr(ell, fs=1.0, rc=0.43): 
    ''' Calculate Del P^uncorr(k) for k values specified by Nseries P(k). 
    Write [k, Del P^uncorr(k)] to pickle file.
    '''
    # k values from Nseries P(k) 
    true_pk_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/', 
        'POWER_Q_CutskyN1.fidcosmo.dat.grid960.P020000.box3600'])
    k_true = np.loadtxt(true_pk_file, unpack = True, usecols =[0])
        
    uncorrdelPk = fourier_corr.delP_uncorr(
            k_true,             # k 
            ell,                # ell 
            fs=fs,              # fs
            rc=rc               # rc
            )
    pickle_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/',
        'uncorrdelP', str(ell), 'k', 
        '.fs', str(round(fs,2)), 
        '.rc', str(round(rc,2)), 
        '.p'
        ])
    pickle.dump([k_true, uncorrdelPk], open(pickle_file, 'wb'))
    return None

def delPcorr_0_qmax(mock, qmax, fs=1.0, rc=0.43):
    ''' Calculate Del P^corr integrated from q = 0 to q = q_max. 
    Dump the [k, Del P^corr(k)] to pickle file 
    '''
    # calculate average true power spectrum
    if mock == 'nseries':       # Nseries
        n_mock = 84
        avgpk_name = 'AVG_POWER_Q_CutskyN.fidcosmo.dat.grid960.P020000.box3600'
        data_dir = '/mount/riachuelo1/hahn/power/Nseries/'
        for i_mock in xrange(1, n_mock+1):
            true_pk_file = ''.join([
                data_dir, 
                'POWER_Q_CutskyN', str(i_mock), '.fidcosmo.dat.grid960.P020000.box3600'])

            tr_k, tr_p0k_i, tr_p2k_i, tr_p4k_i = np.loadtxt(
                        true_pk_file, 
                        unpack = True, 
                        usecols =[0,1,2,3] 
                        )
            if i_mock == 1: 
                tr_p0k = tr_p0k_i
                tr_p2k = tr_p2k_i
                tr_p4k = tr_p4k_i
            else: 
                tr_p0k += tr_p0k_i
                tr_p2k += tr_p2k_i
                tr_p4k += tr_p4k_i 
        tr_specs = [spec/np.float(n_mock) for spec in tr_specs]

    elif mock == 'nseriesbox':  # Nseries Box
        n_mock = 7
        avgpk_name = 'AVG_power3600z_BoxN.dat'
        data_dir = '/mount/riachuelo1/hahn/power/Nseries/Box/'
        for i_mock in xrange(1, n_mock+1):
            true_pk_file = ''.join([data_dir, 'power3600z_BoxN', str(i_mock), '.dat'])

            tr_k, tr_p0k_i, tr_p2k_i, tr_p4k_i, tr_p6k_i, tr_p8k_i, tr_p10k_i = np.loadtxt(
                        true_pk_file, 
                        unpack = True, 
                        usecols =[0,-1,2,3,4,5,6] 
                        )

            if i_mock == 1: 
                tr_p0k = tr_p0k_i
                tr_p2k = tr_p2k_i
                tr_p4k = tr_p4k_i
                tr_p6k = tr_p6k_i
                tr_p8k = tr_p8k_i
                tr_p10k = tr_p10k_i
            else: 
                tr_p0k += tr_p0k_i
                tr_p2k += tr_p2k_i
                tr_p4k += tr_p4k_i
                tr_p6k += tr_p6k_i
                tr_p8k += tr_p8k_i
                tr_p10k += tr_p10k_i

        notnorm_tr_specs = [tr_p0k, tr_p2k, tr_p4k, tr_p6k, tr_p8k, tr_p10k]
        tr_specs =  [(2.0*np.pi)**3/np.float(n_mock) * tr_spec for tr_spec in notnorm_tr_specs]
    else: 
        raise NotImplemented

    if qmax > tr_k[-1]: 
        raise ValueError("qmax can not go beyond max(k)")
    
    for ell in [0, 2]: # calculate the correlated portion of delP_ell(k)
        lps, corrdelPk_lp = fourier_corr.delP_corr_qmax(
                tr_k,               # k 
                tr_specs,           # [p0k, p2k, p4k, ..] 
                ell,                # ell
                q_max=qmax, 
                fs=fs,              # fs 
                rc=rc,              # rc 
                lp_comp=True
                )
    
        corrdelPk = np.zeros(len(tr_k))
        for i_lp, lp in enumerate(lps):     # save l' components to pickle file 

            corrdelPk_pickle_file = ''.join([
                data_dir, 'corrdelP', str(ell), 'k_lp', str(lp), 
                '_qmax', 
                str(round(qmax,2)), '_', avgpk_name, 
                '.rc', str(round(rc, 2)), 
                '.fs', str(round(fs, 2)), 
                '.p'
                ])
            pickle.dump([tr_k, corrdelPk_lp[i_lp]], open(corrdelPk_pickle_file, 'wb'))
            corrdelPk += corrdelPk_lp[i_lp]
        # save Sum_l' to pickle file 
        corrdelPk_pickle_file = ''.join([
            data_dir, 
            'corrdelP', str(ell), 'k', 
            '_qmax', str(round(qmax,2)),'_',avgpk_name, 
            '.rc', str(round(rc, 2)), 
            '.fs', str(round(fs, 2)), 
            '.p'
            ])
        pickle.dump([tr_k, corrdelPk], open(corrdelPk_pickle_file, 'wb'))
    
    msg = 'q_max = '+str(qmax)+'\n'+'mock = '+mock
    codif.notif(subject='delPcorr_0_qmax script finished', message=msg)
    return None

def k_n_coeff_DelPcorr_ktrust_qmax(ktrust, ell=2, fs=0.6, rc=0.43, Nthreads=1, 
        fold=10, rebin=20, order=10):
    ''' Calculate the coefficient of the k^n term of the polynomial 
    expansion of delP_corr integrated from ktrust to q_max. 
    In other words

        \Del P^{corr} |_{q = k_{trust}}^{q = q_{max}}

    Only coded for NseriesBox.

    Parameter
    ---------
    ktrust : float
        The redshift up till we reliably trust the perturbation theory model 
    ell : int
        l = 0 for monopole and l = 2 for quadrupole
    rc : float
        rc is the fiber collision radius. At z~0.35, 0.55, 0.7, rc = 0.35, 0.55, 0.7. 
    fs : float
        fs amplitude of the tophat.

    Notes
    -----
    - shot noise has to be corrected for Nseries Box
    '''
    # list of interpolation function for each P_ell'(k)
    alpha = -fs * rc
    n_ks = [2 * nn for nn in range(order/2+1)]
    f_interps, tr_specs = [], []
    for n_k in n_ks: 
        tr_k, avg_plk = average_Plk(n_k, fold=fold, rebin=rebin)
        tr_specs.append(avg_plk)
        f_interps.append(interp1d(tr_k, avg_plk))#, kind='cubic'))
    
    integ_kwargs_list = []   
    for i_spec, tr_spec in enumerate(tr_specs):
        integ_kwargs_list.append({
                    'ell': ell, 
                    'n_k': n_ks[i_spec],
                    'f_interps': f_interps, 
                    'k_min': tr_k[0], 
                    'k_max': tr_k[-1], 
                    'ktrust': ktrust,
                    'rc': rc
                    })
    
    if Nthreads > 1:   # multiprocessing
        pool = Pewl(processes=Nthreads)
        mapfn = pool.map
        arglist = [[n_ks[ii], integ_kwargs_list[ii]] for ii in range(len(tr_specs))]
        
        n_k_coeff = mapfn(k_n_coeff_DelPcorr_integral, [arg for arg in arglist])

        pool.close()
        pool.terminate()
    else: 
        n_k_coeff = []
        for i_n_k, n_k in enumerate(n_ks): 
            #print 'k^', n_k, ' term' 
            coeff = k_n_coeff_DelPcorr_integral([n_k, integ_kwargs_list[i_n_k]])
            n_k_coeff.append(coeff)
    
    data_list = [np.array(n_ks), alpha * np.array(n_k_coeff)]
    
    fold_str = ''
    if fold: 
        fold_str = ''.join([
            '.', str(fold), 'fold', 
            '.', str(rebin), 'rebin'
            ])

    coeff_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'k^n_coeff.corrDelP_ktrust_qmax',
        '.ell', str(ell), 
        '.ktrust', str(round(ktrust,2)), 
        '.rc', str(round(rc,2)), 
        '.fs', str(round(fs,2)),
        fold_str, 
        '.order', str(order), 
        '.dat'])
    print coeff_file
    np.savetxt(coeff_file, 
            (np.vstack(np.array(data_list))).T, 
            fmt=['%10.5f', '%10.5f'], 
            delimiter='\t', 
            header='Columns : n_k, coeff'
            )
    
    #msg = 'ell = '+str(ell)+'\nk_trust = '+str(ktrust)
    #codif.notif(subject='delPcorr_ktrust_qmax k^n coeff script finished', message=msg)
    return None

def k_n_coeff_DelPcorr_integral(arg): 
    n_k = arg[0] 
    integ_kwargs = arg[1]
    
    rc = integ_kwargs['rc']
    ell = integ_kwargs['ell']
    f_interps = integ_kwargs['f_interps']
    k_min = integ_kwargs['k_min']
    k_max = integ_kwargs['k_max']
    ktrust = integ_kwargs['ktrust']

    # calculate the coefficients 
    coeff_integrand = lambda qq: \
            (J1(qq * rc)/(qq**n_k)) * \
            np.sum([
                poly_exp_coeff(ell, n_k)[ii] * \
                        pq_noextrap(qq, f_interps[ii], q_min=k_min, q_max=k_max)
                for ii in range(len(f_interps))
                ])
    
    coeff_int, coeff_int_err = quad_int(coeff_integrand, ktrust, k_max, epsrel=0.001)
    #options = {'limit': 200}
    #coeff_int, coeff_int_err = nquad(coeff_integrand, [[ktrust, k_max]], opts=[options])

    return coeff_int

def poly_exp_coeff(ell, n_k):
    ''' coefficients for the polynomial expansion. 
    '''
    if ell == 0: 
        if n_k == 0: 
            Clp = [1., -0.5, 0.5*(3./4.), 0.5*(-5./8.), 0.5*(35./64.), 0.5*(-63./128.), 
                    0.5*(231./512.), 0.5*(-429./1024.), 
                    0.5*(6435./16384.), # 16
                    0.5*(-12155./32768.), # 18
                    0.5*(46189./131072.), # 20
                    0.5*(-88179./262144.), # 22
                    0.5*(676039./2097152.)  # 24
                    ]
        elif n_k == 2: 
            Clp = [0., 0.5, 0.5*(-5./2.), 0.5*(35./8.), 0.5*(-105./16.), 0.5*(1155./128.), 
                    0.5*(-6006./512.), 0.5*(15015./1024.), 
                    0.5*(-291720./16384.), # 16
                    0.5*(692835./32768.),  # 18
                    0.5*(-3233230./131072.),  # 20
                    0.5*(7436429./262144.),  # 22
                    0.5*(-67603900./2097152.)    # 24
                    ]
        elif n_k == 4: 
            Clp = [0., 0., 0.5*(7./4.), 0.5*(-63./8.), 0.5*(693./32.), 0.5*(-3003./64.),
                    0.5*(45045./512.), 0.5*(-153153./1024.), 
                    0.5*(3879876./16384.),      # 16
                    0.5*(-11639628./32768.),    # 18
                    0.5*(66927861./131072.),   # 20
                    0.5*(-185910725./262144.),   # 22
                    0.5*(2007835830./2097152.)    # 24
                    ]
        elif n_k == 6: 
            Clp = [0., 0., 0., 0.5*(33./8.), 0.5*(-429./16.), 0.5*(6435./64.), 
                    0.5*(-145860./512.), 0.5*(692835./1024.),
                    0.5*(-23279256./16384.), # 16
                    0.5*(89237148./32768.),   # 18
                    0.5*(-637408200./131072.),   # 20
                    0.5*(2151252675./262144.),   # 22
                    0.5*(-27727256700./2097152.)    # 24
                    ]
        elif n_k == 8: 
            Clp = [0., 0., 0., 0., 0.5*(715./64.), 0.5*(-12155./128.), 
                    0.5*(230945./512.), 0.5*(-1616615./1024.),
                    0.5*(74364290./16384.), # 16
                    0.5*(-371821450./32768.), # 18
                    0.5*(3346393050./131072.), # 20
                    0.5*(-13863628350./262144.), # 22
                    0.5*(214886239425./2097152.)]   # 24
        elif n_k == 10: 
            Clp = [0., 0., 0., 0., 0., 0.5*(4199./128.),
                    0.5*(-176358./512.), 0.5*(2028117./1024.),
                    0.5*(-135207800./16384.),   # 16
                    0.5*(912652650./32768.),    # 18, 
                    0.5*(-10586770740./131072.), #20
                    0.5*(54698315490./262144.), # 22
                    0.5*(-1031453949240./2097152.) # 24
                    ]
        elif n_k == 12: 
            Clp = [0., 0., 0., 0., 0., 0.,
                    0.5*(52003./512.), 0.5*(-1300075./1024.),
                    0.5*(140408100./16384.),    # 16
                    0.5*(-1357278300./32768.),  # 18 
                    0.5*(21037813650./131072.), # 20
                    0.5*(-138849570090./262144.), # 22
                    0.5*(3239823302100./2097152.) # 24
                    ]
        elif n_k == 14: 
            Clp = [0., 0., 0., 0., 0., 0.,
                    0., 0.5*(334305./1024.), 
                    0.5*(-77558760./16384.),    # 16
                    0.5*(1202160780./32768.),   # 18
                    0.5*(-26447537160./131072.), # 20
                    0.5*(231415950150./262144.), # 22
                    0.5*(-6849912124440./2097152.) # 24
                    ]
        elif n_k == 16: 
            Clp = [0., 0., 0., 0., 0., 0.,
                    0., 0., 
                    0.5*(17678835./16384.),     # 16
                    0.5*(-583401555./32768.),   # 18
                    0.5*(20419054425./131072.), # 20
                    0.5*(-251835004575./262144.), # 22
                    0.5*(9821565178425./2097152.) #24
                    ]
        elif n_k == 18:
            Clp = [0., 0., 0., 0., 0., 0.,
                    0., 0., 
                    0., # 16
                    0.5*(119409675./32768.),     # 18
                    0.5*(-8836315950./131072.), # 20 
                    0.5*(172308161025./262144.), # 22
                    0.5*(-9419512802700./2097152.) # 24
                    ]
        elif n_k == 20: 
            Clp = [0., 0., 0., 0., 0., 0.,
                    0., 0., 
                    0., # 16
                    0., # 18
                    0.5*(1641030105./131072.), # 20
                    0.5*(-67282234305./262144.), # 22
                    0.5*(5786272150230./2097152.) #24
                    ]
        elif n_k == 22: 
            Clp = [0., 0., 0., 0., 0., 0.,
                    0., 0., 
                    0., 0., 0.,
                    0.5*(11435320455./262144.), # 22
                    0.5*(-2058357681900./2097152.) # 24
                    ]
        elif n_k == 24: 
            Clp = [0., 0., 0., 0., 0., 0.,
                    0., 0., 
                    0., 0., 0., 0., 
                    0.5*(322476036831./2097152.) # 24
                    ]
        else: 
            raise NotImplementedError
    elif ell == 2: 
        if n_k == 0: 
            Clp = [0., 0., 0., 0., 0., 0.,
                    0., 0., 
                    0., # 16
                    0., # 18
                    0., # 20 
                    0., # 22 
                    0. # 24 
                    ]
        elif n_k == 2: 
            Clp = [0., 1., -2.5, 2.5*7./4., 2.5*(-21./8.), 2.5*(231./64.), 
                    2.5*(-3003./640.), 2.5*(3003./512.),
                    2.5*(-7293./1024.), # 16
                    2.5*(138567./16384.), # 18
                    2.5*(-323323./32768.), # 20 
                    2.5*(7436429./655360.), # 22 
                    2.5*(-3380195./262144.) # 24 
                    ]
        elif n_k == 4: 
            Clp = [0., 0., 2.5, 2.5*(-9./2.), 2.5*(99./8.), 2.5*(-429./16.), 
                    2.5*(32175./640.), 2.5*(-43758./512.), 
                    2.5*(138567./1024.), # 16
                    2.5*(-3325608./16384.), # 18
                    2.5*(9561123./32768.), # 20 
                    2.5*(-265586750./655360.), # 22 
                    2.5*(143416845./262144.) # 24 
                    ]
        elif n_k == 6: 
            Clp = [0., 0., 0., 2.5*(11./4.), 2.5*(-143./8.), 2.5*(2145./32.), 
                    2.5*(-121550./640.), 2.5*(230945./512.), 
                    2.5*(-969969./1024.), # 16
                    2.5*(29745716./16384.), # 18
                    2.5*(-106234700./32768.), # 20 
                    2.5*(3585421125./655360.), # 22 
                    2.5*(-2310604725./262144.) # 24 
                    ]
        elif n_k == 8: 
            Clp = [0., 0., 0., 0., 2.5*(65./8.), 2.5*(-1105./16.), 
                    2.5*(209950./640.), 2.5*(-587860./512.), 
                    2.5*(3380195./1024.), # 16
                    2.5*(-135207800./16384.), # 18
                    2.5*(608435100./32768.), # 20 
                    2.5*(-25206597000./655360.), # 22 
                    2.5*(19535112675./262144.) # 24 
                    ]
        elif n_k == 10: 
            Clp = [0., 0., 0., 0., 0., 2.5*(1615./64.),
                    2.5*(-169575./640.), 2.5*(780045./512.),
                    2.5*(-6500375./1024.), # 16
                    2.5*(351020250./16384.), # 18
                    2.5*(-2035917450./32768.), # 20 
                    2.5*(105189068250./655360.), # 22 
                    2.5*(-99178264350./262144.) # 24 
                    ]
        elif n_k == 12: 
            Clp = [0., 0., 0., 0., 0., 0.,
                    2.5*(52003./640.), 2.5*(-520030./512.), 
                    2.5*(7020405./1024.), # 16
                    2.5*(-542911320./16384.), # 18
                    2.5*(4207562730./32768.), # 20 
                    2.5*(-277699140180./655360.), # 22 
                    2.5*(323982330210./262144.) # 24 
                    ]
        elif n_k == 14:
            Clp = [0., 0., 0., 0., 0., 0., 
                    0., 2.5*(137655./512.),
                    2.5*(-3991995./1024.), # 16
                    2.5*(495007380./16384.), # 18
                    2.5*(-5445081180./32768.), # 20 
                    2.5*(476444603250./655360.), # 22 
                    2.5*(-705138012810./262144.) # 24 
                    ]
        elif n_k == 16: 
            Clp = [0., 0., 0., 0., 0., 0., 
                    0., 0., 
                    2.5*(930465./1024.), # 16
                    2.5*(-245642760./16384.), # 18
                    2.5*(4298748300./32768.), # 20 
                    2.5*(-530178957000./655360.), # 22 
                    2.5*(1033848966150./262144.) # 24 
                    ]
        elif n_k == 18: 
            Clp = [0., 0., 0., 0., 0., 0., 
                    0., 0., 0., 
                    2.5*(51175575./16384.), # 18
                    2.5*(-1893496275./32768.), # 20 
                    2.5*(369231773625./655360.), # 22 
                    2.5*(-1009233514575./262144.) # 24 
                    ]
        elif n_k == 20: 
            Clp = [0., 0., 0., 0., 0., 0., 
                    0., 0., 0., 0.,
                    2.5*(356745675./32768.), # 20 
                    2.5*(-146265726750./655360.), # 22 
                    2.5*(628942625025./262144.) # 24 
                    ]
        elif n_k == 22: 
            Clp = [0., 0., 0., 0., 0., 0., 
                    0., 0., 0., 0., 0.,
                    2.5*(25157705001./655360.), # 22 
                    2.5*(-226419345009./262144.) # 24 
                    ]
        elif n_k == 24: 
            Clp = [0., 0., 0., 0., 0., 0., 
                    0., 0., 0., 0., 0., 0.,
                    2.5*(35830670759./262144.) # 24 
                    ]
        else: 
            raise NotImplementedError
    else:
        raise NotImplementedError('Only quadrupole correction implemented') 

    return Clp

def pq_noextrap(q, f_interp, q_min=0.0003, q_max=4.34): 
    if (q > q_min) and (q <= q_max): 
        return f_interp(q)
    elif (q <= q_min): 
        return 0.0
    elif (q > q_max): 
        return 0.0

def average_Plk(ell, fold=10, rebin=20): 
    ''' Average power spectrum multipole 
    '''
    n_mock = 3 
    for i in range(1, n_mock+1): 
        k, plk = nseriesbox_Plk(ell, i, Ngrid=3600, fold=fold, rebin=rebin)
        
        if i == 1: 
            avg_k = k 
            avg_plk = plk
        else: 
            avg_k += k 
            avg_plk += plk
    
    avg_k /= np.float(n_mock)
    avg_plk *= (2.0*np.pi)**3/np.float(n_mock)
    
    return [avg_k, avg_plk]

def nseriesbox_Plk(ell, imock, Ngrid=3600, fold=10, rebin=20):
    ''' Nseries Box Power spectrum multipoles P_l(k)
    '''
    i_col = ell/2+1
    if ell == 0: 
        i_col = -1 

    pk_file = ''.join([
        '/mount/riachuelo1/hahn/power/Nseries/Box/', 
        'plk3600z_BoxN', str(imock), '.dat'
        ])

    k, plk, co = np.loadtxt(pk_file, unpack=True, usecols=[0, i_col,-2])
    k_max = k[-1]

    if fold is not None:
        # '/export/masa1/rs123/BOSS/', 
        # 'power3600z_BoxN', str(imock), '_', str(fold), 'xfold.dat'
        fold_file = ''.join([
            '/mount/riachuelo1/hahn/power/Nseries/Box/', 
            'plk3600z_BoxN', str(imock), '.10xfold.dat' 
            ])
        k_f, plk_f, co_f = np.loadtxt(fold_file, unpack=True, usecols=[0, i_col,-2])
        plk_f *= fold**3
   
        fold_range = np.where(k_f > k_max)

        k_f = k_f[fold_range]
        plk_f = plk_f[fold_range]
        co_f = co_f[fold_range]

        k_rebin = np.logspace(np.log10(k_max), np.log10(k_f[-1]), rebin)
        avg_k_f, avg_plk_f = [], []
        for i_k in range(1, len(k_rebin)-1): 
            k_logbin_lo = k_rebin[i_k-1]
            k_logbin_hi = k_rebin[i_k+1]

            within = np.where((k_f >= k_logbin_lo) & (k_f < k_logbin_hi))
            if len(within[0]) == 0: 
                continue 

            avg_k_f.append(np.sum(k_f[within] * co_f[within])/np.sum(co_f[within]))
            avg_plk_f.append(np.sum(plk_f[within] * co_f[within])/np.sum(co_f[within]))

        avg_k_f = np.array(avg_k_f)
        avg_plk_f = np.array(avg_plk_f)
        
        k_tot = np.concatenate([k, avg_k_f]).T
        plk_tot = np.concatenate([plk, avg_plk_f]).T

    else: 
        k_tot = k
        plk_tot = plk

    return [k_tot, plk_tot]



if __name__=='__main__': 
    for n_k in [6, 8, 10, 12, 14, 16, 18, 20, 22, 24]: 
        k_n_coeff_DelPcorr_ktrust_qmax(
                0.3, ell=0, fs=0.6, rc=0.43, fold=10, rebin=50, Nthreads=1, order=n_k)
        k_n_coeff_DelPcorr_ktrust_qmax(
                0.3, ell=2, fs=0.6, rc=0.43, fold=10, rebin=50, Nthreads=1, order=n_k)
