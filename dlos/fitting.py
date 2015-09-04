"""

Fitting the peak of the line-of-sight displacement distribution


"""
from dlos import Dlos

def dlospeak_fit(
        dlos, 
        fit = 'gauss', 
        peak_range = [-15.0, 15.0]
        ): 
    """ Fit the peak of the dLOS distribution to specified 
    fit type. Function uses MPfit. 
    """
    
    # dummy catalog correction dictionary
    dum_cat_corr = {
            'catalog': {'name': 'nseries', 'n_mock': 1}, 
            'correction': {'name': 'upweight'}
            }

    dlosclass = Dlos(dum_cat_corr)

    # Freedman Diaconis binsize
    fd_binsize = dlosclass.fd_binsize(
            dlos = dlos, 
            peak_range = peak_range
            )
    
    dlos_mid, dlos_dist = dlosclass.dlos_dist(
            binsize = fd_binsize
            )
    
    # only fit the peak of the dLOS distribution  
    inpeak = np.where(
            (dlos_mid < peak_range[1]) & 
            (dlos_mid > peak_range[0])
            ) 

    amp_guess = np.mean( dlos_dist[ np.where( (dlos_mid < 2.0) & (dlos_mid > -2.0) ) ] )

    p0 = [amp_guess, 5.0]    # initial guess

    fa = {'x': dlos_mid[inpeak], 'y': dlos_dist[inpeak]}
    if fit == 'expon': 
        mpfit_func = mpfit_peak_expon
    elif fit == 'gauss': 
        mpfit_func = mpfit_peak_gauss
    else: 
        raise NotImplementedError("Fit not yet coded") 

    fit_param = mpfit.mpfit(mpfit_func, p0, functkw=fa)
        
    best_amp = fit_param.params[0]
    best_sigma = fit_param.params[1]

    if fit == 'gauss': 
        bestfit_func_integral = np.sqrt(2.0 * np.pi) * best_amp * best_sigma 
    else: 
        raise NotImplementedError()
        
    # compute fpeak 
    fpeak_xmin = -3.0*sigma
    fpeak_xmax = 3.0*sigma
    fpeak_xrange = (xmid >= fpeak_xmin) & (xmid < fpeak_xmax) # +/- 3 sigma 

    if fit.lower() == 'expon': 
        fpeak = (np.sum(peak_expon(xmid[fpeak_xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    elif fit.lower() == 'gauss': 
        fpeak = (np.sum(peak_gauss(xmid[fpeak_xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    else: 
        raise NameError("Fit not yet coded") 
        
    fpeak_dist = np.float( len(dlos[ (dlos > fpeak_xmin) & (dlos < fpeak_xmax) ]) )\
            /np.float(len(dlos))  # computed from the distribution

    #print 'fpeak from fitting = ', fpeak
    #print 'fpeak from distrib = ', fpeak_dist 

    fit_dict = {'amp': bestfit_amp, 'sigma': sigma, 'fpeak': fpeak_dist} 

    return dlos_hist, xmid, fit_dict

def combined_dlos_fit(n_mocks, fit='gauss', sanitycheck=False, clobber=False, **cat_corr):
    ''' Fits combined dLOS distribution to specified function type 
    using Freedman-Diaconis binsizes and MPfit

    Parameters
    ----------
    n_mocks : number of mock dLOS files to include 
    fit : fitting function ('gauss', 'expon', etc) 
    sanitycheck : testing flag to plot the dLOS distribution and fit
    cat_corr : catalog and correction dictionary 
    clobber : re-make dLOS files (True or False)  
    
    Returns 
    ------- 
    sigma : sigma value of gaussian/exponential function  
    fpeak : peak fraction of the dLOS distribution 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
   
    # Read dLOS data 
    if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz'): # LasDamasGeo ----------------
        for i_mock in range(1, n_mocks+1): 
            for letter in ['a', 'b', 'c', 'd']: 
                # individual catalog_correction dictonary 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock 
                i_catalog['letter'] = letter 
                i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
                
                los_disp_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
                # Combine dLOS values 
                try: 
                    combined_dlos
                except NameError: 
                    combined_dlos = los_disp_i.dlos
                else: 
                    combined_dlos = np.concatenate([combined_dlos, los_disp_i.dlos]) 

    elif catalog['name'].lower() in ('tilingmock'): # Tiling Mock --------------------------

        # import DLOS values for mock 
        los_disp = dlos(**cat_corr)

        combined_dlos = los_disp.dlos
        n_mocks = 1

    elif 'bigmd' in catalog['name'].lower():        # BigMD --------------------------------
        # import DLOS values for mock 
        los_disp = dlos(**cat_corr)

        combined_dlos = los_disp.dlos
        #rand_sub = np.random.randint(len(los_disp.dlos), size=2500)
        #combined_dlos = (los_disp.dlos)[rand_sub]
        n_mocks = 1

    elif catalog['name'].lower() in ('qpm', 'patchy', 'nseries'):   # QPM and PATCHY 
        for i_mock in range(1, n_mocks+1): 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            los_disp_i = dlos(clobber=clobber, **i_cat_corr)   # import DLOS of each mock 

            # Combine dLOS values 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i.dlos
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i.dlos]) 
    
    elif 'cmass' in catalog['name'].lower():        # CMASS ---------------------------------
        # import DLOS values for mock 
        los_disp = dlos(**cat_corr)

        combined_dlos = los_disp.dlos
        n_mocks = 1
    else: 
        raise NameError("not yet coded") 

    # Determine appropriate binsize for dist (not too important)
    # Create histogram for combined dLOS values  (binsize is just guessed)
    x_min = -1000.0
    x_max = 1000.0
    binsize = 0.1 
    n_bins = int((x_max-x_min)/binsize) 
    
    guess_dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
    guess_xlow = mpc_binedges[:-1]
    guess_xhigh = mpc_binedges[1:] 
    guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
    
    # determine appropriate bin size using Freedman-Diaconis Rule 
    peak_range = (guess_xmid >= 0.0) & (guess_xmid < 15.0) 
    dlos_cumu = (guess_dlos_hist[peak_range]).cumsum() 
    n_sample = dlos_cumu[-1]

    iqr_index = fc_util.find_nearest(dlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
    
    binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
    print 'Freedman-Diaconis binsize = ', binsize
    #-------------------------------------------------------------------------------------
    # recompute histogram using freedman-diaconis binsize 
    n_bins = int((x_max-x_min)/binsize) 

    dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])

    # fitting ------------------------------------------------------------------------------
    #if fit.lower() == 'gauss': 
    #elif fit.lower() == 'expon': 
    #    peak_xrange = (xmid > -25.0) & (xmid < 25.0) 
    #else: 
    #    raise NameError('Error')
    peak_xrange = (xmid > -25.0) & (xmid < 25.0)            # rough peak range 

    # amplitude of dLOS distribution 
    dlos_amp = np.mean(dlos_hist[(xmid >= -1.0) & (xmid < 1.0)])
    #print 'Amplitude = ', dlos_amp, np.max(dlos_hist)
  
    # MPFIT ------------------------------------------------------------------------------------
    p0 = [dlos_amp, 5.0]            # initial guess

    fa = {'x': xmid[peak_xrange], 'y': dlos_hist[peak_xrange]}
    if fit.lower() == 'expon': 
        peak_pars = mpfit.mpfit(mpfit_peak_expon, p0, functkw=fa, nprint=0)
    elif fit.lower() == 'gauss': 
        peak_pars = mpfit.mpfit(mpfit_peak_gauss, p0, functkw=fa, nprint=0)
    else: 
        raise NameError("Fit not yet coded") 
    
    bestfit_amp = peak_pars.params[0]
    sigma = peak_pars.params[1]
    print fit.lower(), ' Best FIt Sigma = ', sigma
    
    #--------------------------------------------------------------------------------------------------------------------
    # compute fpeak 
    fpeak_xmin = -3.0*sigma
    fpeak_xmax = 3.0*sigma
    xrange = (xmid >= fpeak_xmin) & (xmid < fpeak_xmax) 
    #fpeak = (np.sum(peak(xmid[xrange], popt[0])))/np.float(np.sum(dlos_hist)) 
    if fit.lower() == 'expon': 
        fpeak = (np.sum(peak_expon(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    elif fit.lower() == 'gauss': 
        fpeak = (np.sum(peak_gauss(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    else: 
        raise NameError("Fit not yet coded") 
    
    fpeak_dist = np.float(len(combined_dlos[(combined_dlos > fpeak_xmin) & (combined_dlos < fpeak_xmax)]))/np.float(len(combined_dlos))  # computed from the distribution

    print 'fpeak from fitting = ', fpeak
    print 'fpeak from distrib = ', fpeak_dist 
    
    #fpeak_dist = np.float(len(combined_dlos[(combined_dlos > -30.0) & (combined_dlos < 30.0)]))/np.float(len(combined_dlos))  # computed from the distribution

    #print 'fpeak from distrib = ', fpeak_dist 

    if sanitycheck == True:             # plot the distributions and fits for santicy check 
        prettyplot() 
        pretty_colors = prettycolors() 
        fig = plt.figure(1, figsize=[14,5]) 
        sub = fig.add_subplot(111) 

        # fitting labels 
        if fit.lower() == 'expon': 
            fit_label = "Exponential "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
        else: 
            fit_label = "Gaussian "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
    
        #guess_scale = binsize/0.1
        #sub.plot(guess_xmid, guess_scale*guess_dlos_hist, lw=2, color=pretty_colors[-1], label=r"$d_{LOS}$ binsize $=0.1$") 

        sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[0], label=r"$d_{LOS}$ binsize $="+str(binsize)+"$")         # plot DLOS distribution

        # plot best fit 
        if fit.lower() == 'expon': 
            sub.plot(xmid, peak_expon(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
        elif fit.lower() == 'gauss': 
            sub.plot(xmid, peak_gauss(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
        else: 
            raise NameError("Fit not yet coded") 

        # indicate the fpeak range (for sanity check purposes)
        sub.vlines(fpeak_xmin, 0.0, 10.0*np.max(dlos_hist), color='k', lw=4)
        sub.vlines(fpeak_xmax, 0.0, 10.0*np.max(dlos_hist), color='k', lw=4)

        sub.text(-1.0*sigma, 0.25*np.max(dlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
        sub.text(-1.0*sigma, 0.2*np.max(dlos_hist), r"$f_{peak,dist} = "+str(fpeak_dist)+"$") 
        sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
        sub.set_xlim([-50.0, 50.0])
        #sub.set_xlim([-5.0, 5.0])
        sub.set_ylim([0.0, 1.25*np.max(dlos_hist)])
        sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 
    
        if correction['name'].lower() == 'bigfc': 
            bigfc_flag = '_bigfc'
        else: 
            bigfc_flag = ''

        fig_dir = 'figure/'
        fig_file = ''.join([fig_dir, 
            catalog['name'].lower(), '_', str(n_mocks), 
            'mocks_combined_dlos_peakfit_', fit.lower(), bigfc_flag, '.png'])
        fig.savefig(fig_file, bbox_inches="tight")
        fig.clear() 

    return [sigma, fpeak]

#---- Fitting -----
def fit_linear(x, p): 
    ''' Linear function y = a * x + b 
    p[0] = a , p[1] = b
    '''
    return (p[0] * x) + p[1]

def peak_expon(x, p): 
    ''' Exponential function for the peak 
    '''
    return p[0]*np.exp(-np.abs(x)/p[1])

def peak_gauss(x, p): 
    """ Exponential function for the peak 
    """
    return p[0]*np.exp(-0.5*x**2/(p[1])**2)

# --- MPfit ---
def mpfit_linear(p, fjac=None, x=None, y=None, err=None): 
    model = fit_linear(x, p) 
    status = 0 
    if err == None: 
        err = np.array([1.0 for i in range(len(x))])
    return([status, (y-model)/err]) 

def mpfit_peak_expon(p, fjac=None, x=None, y=None): 
    model = peak_expon(x, p) 
    status = 0 
    return([status, (y-model)]) 

def mpfit_peak_gauss(p, fjac=None, x=None, y=None): 
    model = peak_gauss(x, p) 
    status = 0 
    return([status, (y-model)]) 
