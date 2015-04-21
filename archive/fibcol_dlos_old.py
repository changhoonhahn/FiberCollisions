# ARCHIVE OF OLD STUFF
"""
def build_dlos_corr_cosmo(**cat_corr): 
    '''
    build DLOS using correct cosmology
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    if catalog['name'].lower() != 'qpm': raise NameError('asdfasdfasdf') 

    mock_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
        'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz']) 
    print mock_file
    los_disp = dlos(readdata=False, **cat_corr)  
    dlos_file = '.corr.cosmo.'.join((los_disp.file_name).rsplit('.', 1)) 
    
    build_dlos_idl_cmd = ''. join(['idl -e ', '"', 'build_fibcoll_dlos_cosmo, ', "'", catalog['name'], "','", mock_file, "','", dlos_file, "'", '"'])
    print build_dlos_idl_cmd 
    os.system(build_dlos_idl_cmd) 
"""

"""
def combine_dlos_correct_cosmology(fit, **cat_corr): 
    '''
    combined dlos distribution fitting using *correct cosmology*
    '''
    catalog = cat_corr['catalog'] 
    correction = {'name': 'true'} 
    
    n_mocks = 0 
    if catalog['name'].lower() == 'qpm': 
        for i_mock in range(1, 46):             # hardcoded for only 10 files 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            los_disp_rsd_i = dlos(readdata=False, **i_cat_corr)         # import DLOS values from each mock 
            los_disp_i_file = '.corr.cosmo.'.join((los_disp_rsd_i.file_name).rsplit('.', 1))       # hardcoded to read in correct cosmology dLOS 

            if os.path.isfile(los_disp_i_file) == False: 
                build_dlos_corr_cosmo(**i_cat_corr) 
            
            los_disp_i = np.loadtxt(los_disp_i_file, unpack=True, usecols=[0]) 

            # Combine dLOS values 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i]) 
            n_mocks = n_mocks+1

    # Get Freedman Diaconis binsizes ------------------------------------------------------------
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
    # -----------------------------------------------------------------------------------------------

    # recompute histogram using appropriate binsize 
    n_bins = int((x_max-x_min)/binsize) 
    dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
    print 'binsize ', binsize
    print n_bins
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
    # amplitude of dLOS distribution 
    dlos_amp = np.mean(dlos_hist[(xmid >= -1.0) & (xmid < 1.0)])
    print dlos_amp, np.max(dlos_hist)
   
    if fit.lower() == 'gauss': 
        peak_xrange = (xmid > -25.0) & (xmid < 25.0)            # ~ +/- 5 sigmas if sigma ~5
    elif fit.lower() == 'expon': 
        peak_xrange = (xmid > -25.0) & (xmid < 25.0) 
    else: 
        raise NameError('Error')

    # MPFit fitting ----------------------------------------------------------------------------------------
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
    print 'bestfit sigma = ', sigma
    #----------------------------------------------------------------------------------------------------------------
    # Compute f_peak
    xrange = (xmid >= -5.0*sigma) & (xmid < 5.0*sigma) 
    #fpeak = (np.sum(peak(xmid[xrange], popt[0])))/np.float(np.sum(dlos_hist)) 

    # from the mpfit
    if fit.lower() == 'expon': 
        fpeak = (np.sum(peak_expon(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    elif fit.lower() == 'gauss': 
        fpeak = (np.sum(peak_gauss(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    else: 
        raise NameError("Fit not yet coded") 

    print '+/-5 sigma x-range ', -5.0*sigma, '< d < ', 5.0*sigma
    print 'fpeak, calculated from fit = ', fpeak

    fpeak_dist = np.float(len(combined_dlos[(combined_dlos > -5.0*sigma) & (combined_dlos < 5.0*sigma)]))/np.float(len(combined_dlos))  # computed from the distribution
    print 'fpeak, calculated from dist = ', fpeak_dist 
    
    # plotting ------------------------------------------------------------------------------------------------
    # fitting labels 
    if fit.lower() == 'expon': 
        fit_label = "Exponential "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
    else: 
        fit_label = "Gaussian "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"

    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    guess_scale = binsize/0.1
    sub.plot(guess_xmid, guess_scale*guess_dlos_hist, lw=2, color=pretty_colors[-1], label=r"$d_{LOS}$ binsize $=0.1$") 
    sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[0], label=r"$d_{LOS}$ binsize $="+str(binsize)+"$") 

    if fit.lower() == 'expon': 
        sub.plot(xmid, peak_expon(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
    elif fit.lower() == 'gauss': 
        sub.plot(xmid, peak_gauss(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
    else: 
        raise NameError("Fit not yet coded") 

    sub.text(-1.0*sigma, 0.25*np.max(dlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
    sub.text(-1.0*sigma, 0.23*np.max(dlos_hist), r"$f_{peak,dist} = "+str(fpeak_dist)+"$") 
    sub.set_xlabel(r"$d_{LOS}$", fontsize=20) 
    sub.set_xlim([-20.0, 20.0])
    sub.set_ylim([0.0, 1.25*np.max(dlos_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_dlos_corr_cosmo_'+fit.lower()+'_fit.png', bbox_inches="tight")
    fig.clear() 
"""
'''
def plot_curvefit_dlos(catalog='PTHalo', version='v11p0', n_mock=1, fit='expon'): 
    plot the bestfit exponential/gaussian to dLOS histogram 
    dlos xrange is currently hardcoded
    prettyplot()
    dlos_hist = dlos_histogram()
    dlos_hist.readFile(catalog=catalog, version=version, n_mock=n_mock) 
    dlos_hist.get_hist()
    print dlos_hist.file 

    if fit=='expon': 
        popt, pcov = curve_fit(dlos_hist.expon, np.array(dlos_hist.xmid[10000:11000]), np.array(dlos_hist.hist[10000:11000]))
    elif fit=='gauss':
        popt, pcov = curve_fit(dlos_hist.gauss, np.array(dlos_hist.xmid[9500:10500]), np.array(dlos_hist.hist[9500:10500]))
    print popt

    fig = py.figure(1)
    sub = fig.add_subplot(111)
    sub.plot(dlos_hist.xmid, dlos_hist.hist, linewidth=3, label=r'Histogram of $\mathtt{d_{LOS}}$')
    sub.plot(dlos_hist.xmid[10000:10500], dlos_hist.expon(np.array(dlos_hist.xmid[10000:10500]), popt[0]), \
            'r', linewidth=3, label=r'Exponential bestfit with $\mathtt{\sigma='+str(popt[0])+'}$')
    sub.set_xlim([-50, 50])
    ylimit = [0.0, np.max(dlos_hist.hist)*1.25]
    sub.set_ylim(ylimit)
    sub.set_xlabel('Line-of-Sight Displacement (Mpc)')
    sub.set_ylabel('Number of Galaxies') 
    sub.legend(loc='upper left', fontsize=15)
    for d in [20, 30, 40]:
        RMSfrac = float(len(dlos_hist.dlos[(dlos_hist.dlos < d) & (dlos_hist.dlos > -d)]))/float(len(dlos_hist.dlos))*100.0
        caption = r''+str(np.int(RMSfrac))+"$\%$"
        sub.annotate(caption, (float(d), ylimit[1]/10.0), xycoords='data', xytext=(float(d), ylimit[1]/2.0), textcoords='data', \
            arrowprops=dict(arrowstyle="fancy", facecolor='black', connectionstyle="angle3,angleA=0,angleB=-90"),\
            fontsize=20, horizontalalignment='center', verticalalignment='top')

    if dlos_hist.catalog=='LasDamasGeo': 
        file_id = ''.join([dlos_hist.n_mock])
        print file_id
    elif dlos_hist.catalog=='LasDamasBox': 
        print 'asdfas'
        #file_id = str(dlos_hist.filespec[0])
    elif dlos_hist.catalog == 'PTHalo': 
        file_id = ''.join([dlos_hist.version, '_', str(dlos_hist.n_mock)])
    elif (dlos_hist.catalog).lower() == 'qpm': 
        file_id = ''.join([dlos_hist.version, '_', str(dlos_hist.n_mock)])
    else: 
        file_id = ''
    #sub.text(-40.0, ylimit[1]*0.75, ''.join([dlos_hist.catalog,' ', file_id]), fontsize=15) 

    fig_name = ''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', 
        'dLOS_distribution_bestfit_', dlos_hist.catalog, file_id, '_', fit, 'fit.png'])
    print fig_name
    fig.savefig(fig_name)
    fig.clear()
'''    
