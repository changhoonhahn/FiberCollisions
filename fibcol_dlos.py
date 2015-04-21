'''


Line of Sight Displacement Codes

Author(s): ChangHoon Hahn 


'''


import numpy as np
from scipy.optimize import curve_fit
import sys
import os.path
import subprocess
import cosmolopy as cosmos

# -- Local -- 
import fibcol_data as fc_data
import fibcol_utility as fc_util
import mpfit as mpfit
import galaxy_environment as genv

# ------------------------------------------------------------------------------------
# Set up LOS displacement class
# ------------------------------------------------------------------------------------
class dlos:
    def __init__(self, readdata=True, **cat_corr): 
        '''
        Given catalog_correction info, read dLOS values. If dLOS values do not exist, make them
        '''
        catalog = cat_corr['catalog']
        correction = cat_corr['correction']

        # LasDamas Geo ------------------------------------------------------------------------------------------
        if catalog['name'].lower() == 'lasdamasgeo': 

            file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
            #File name 
            file_name = file_dir+'DLOS_sdssmock_gamma_lrgFull_zm_oriana'+str(catalog['n_mock']+100)[1:3]+catalog['letter']+'_no.rdcz.dat'
            self.file_name = file_name 
            
            if readdata == True: 
                if os.path.isfile(file_name) == False:              # if file does not exist then
                    build_dlos(**cat_corr) 
                
                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])

                self.dlos = readin_data[0]
                self.targ_z = readin_data[1]
                self.neigh_z = readin_data[2]
        
        # Tiling Mock ------------------------------------------------------------------------------------------
        elif catalog['name'].lower() == 'tilingmock': 

            file_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'        # directory
            # File name 
            file_name = file_dir+'DLOS_cmass-boss5003sector-icoll012.dat'
            self.file_name = file_name 
            
            if readdata == True: 
                if os.path.isfile(file_name) == False: 
                    # if file does not exist then
                    build_dlos(**cat_corr) 
                
                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])

                self.dlos = readin_data[0]
                self.targ_z = readin_data[1]
                self.neigh_z = readin_data[2]
        
        # QPM ----------------------------------------------------------------------------------------------------
        elif catalog['name'].lower() == 'qpm': 
            file_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/'        # directory
            # File name 
            file_name = ''.join([file_dir, 'DLOS_a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.vetoed.fibcoll.dat']) 
            self.file_name = file_name 
            
            if readdata == True: 
                if os.path.isfile(file_name) == False: 
                    print 'Constructing ', file_name 
                    # if file does not exist then
                    build_dlos(**cat_corr) 

                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3,4,5,6]) 
                self.dlos = readin_data[0]
                self.targ_ra = readin_data[1]
                self.targ_dec = readin_data[2]
                self.targ_z = readin_data[3]
                self.neigh_ra = readin_data[4] 
                self.neigh_dec = readin_data[5] 
                self.neigh_z = readin_data[6]
        
        # CMASS --------------------------------------------------------------------------------------------------------
        elif catalog['name'].lower() == 'cmass': 
            file_dir = '/mount/riachuelo1/hahn/data/'        # directory
            # File name 
            file_name = ''.join([file_dir, 'DLOS_cmass-dr12v4-N-Reid-weights-zlim.dat']) 
            self.file_name = file_name 
            
            if readdata == True: 
                if os.path.isfile(file_name) == False: 
                    print 'Constructing ', file_name 
                    # if file does not exist then
                    build_dlos(**cat_corr) 

                readin_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2]) 
                self.dlos = readin_data[0]
                self.targ_z = readin_data[1]
                self.neigh_z = readin_data[2]

        else: 
            raise NameError('not coded yet') 

    def get_dlos_hist(self, binsize=0.1, normed=False): 
        ''' 
        Compute the histogram for dlos  class
        '''
        x_min = -1000.0
        x_max = 1000.0
        n_bins = int((x_max-x_min)/binsize) 
        dlos_hist, mpc_binedges = np.histogram(self.dlos, bins=n_bins, range=[x_min, x_max], normed=normed)
        # store histogram info 
        self.dlos_hist = dlos_hist 
        self.xlow = mpc_binedges[:-1]
        self.xhigh = mpc_binedges[1:] 
        self.xmid = np.array([0.5*(self.xlow[i]+self.xhigh[i]) for i in range(len(self.xlow))])

    def expon(self, x, sig):
        return np.max(self.dlos_hist)*np.exp(-x/sig)
    def gauss(self, x, sig):
        return np.max(self.dlos_hist)*np.exp(-0.5*x**2/sig**2)

def build_dlos(**cat_corr): 
    '''
    build DLOS for given catalog_correction ID 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
   
    if catalog['name'].lower() == 'qpm': 
        # idl code only works for QPM upweight
        if correction['name'].lower() != 'upweight': 
            correction['name'] = 'upweight'

    elif catalog['name'].lower() == 'lasdamasgeo': 
        # idl code only fit for running LasDamasGeo True Correction 
        if correction['name'].lower() != 'true': 
            correction['name'] = 'true'

    elif catalog['name'].lower() == 'tilingmock': 
        if correction['name'].lower() != 'upweight': 
            correction['name'] = 'upweight'
    
    elif catalog['name'].lower() == 'cmass':  
        if correction['name'].lower() != 'upweight': 
            correction['name'] = 'upweight'

    else: 
        raise NameError('not yet coded')


    mock = fc_data.galaxy_data('data', readdata=True, **cat_corr)
    mock_file = mock.file_name 
    los_disp = dlos(readdata=False, **cat_corr)  
    dlos_file = los_disp.file_name

    #build_dlos_idl_cmd = ''. join(['idl -e ', '"', 'build_fibcoll_dlos, ', "'", catalog['name'], "','", mock_file, "','", dlos_file, "'", '"'])            # incorrect cosmology!
    build_dlos_idl_cmd = ''. join(['idl -e ', '"', 'build_fibcoll_dlos_cosmo, ', "'", catalog['name'], "','", mock_file, "','", dlos_file, "'", '"'])
    
    print build_dlos_idl_cmd 
    os.system(build_dlos_idl_cmd) 

def build_dlos_zreal(**cat_corr): 
    '''
    build DLOS using REAL redshift for given catalog_correction ID 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    if catalog['name'].lower() != 'qpm': raise NameError('asdfasdfasdf') 

    mock_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
        'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz']) 
    print mock_file
    los_disp = dlos(readdata=False, **cat_corr)  
    dlos_file = '.zreal.'.join((los_disp.file_name).rsplit('.', 1)) 
    
    build_dlos_idl_cmd = ''. join(['idl -e ', '"', 'build_fibcoll_dlos_zreal, ', "'", catalog['name'], "','", mock_file, "','", dlos_file, "'", '"'])
    print build_dlos_idl_cmd 
    os.system(build_dlos_idl_cmd) 

def compute_dlos_py(**cat_corr): 
    '''
    compute dlos using cosmolopy and compare dLOS calculated from IDL  
    NO SIGNIFICANT DIFFERENCE
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    # cosmology *DEPENDS ON MOCK CATALOG!* ---------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo': 
        omega_m = 0.25
    elif catalog['name'].lower() == 'qpm': 
        omega_m = 0.31
    elif catalog['name'].lower() == 'tilingmock': 
        omega_m = 0.274
    else: 
        raise NameError('not yet coded!')
    print 'Omega_m = ', omega_m, 'Omega_L = ', 1.0-omega_m      # assumign flatness
    cosmo = {}
    cosmo['omega_M_0'] = omega_m 
    cosmo['omega_lambda_0'] = 1.0 - omega_m 
    cosmo['h'] = 0.7 
    cosmo = cosmos.distance.set_omega_k_0(cosmo) 

    #---------------------------------------------------------------------------------------
    # import dLOS values 
    n_mocks = 0 
    for i_mock in range(1, 11):             # hardcoded to look at only 10 mock DLOS values 
        catalog['n_mock'] = i_mock
        i_cat_corr['catalog'] = catalog
        i_cat_corr['correction'] = correction
        
        los_disp = dlos(readdata=True, **i_cat_corr)                 # dLOS file 

        #read dLOS_IDL, targ_z, neigh_z
        #dLOS_data = np.loadtxt(los_disp.file_name, unpack=True, usecols=[0,1,2])    
        
        targ_Dc = cosmos.distance.comoving_distance(los_disp.targ_z, **cosmo)*cosmo['h']
        neigh_Dc = cosmos.distance.comoving_distance(los_disp.neigh_z, **cosmo)*cosmo['h']

        dLOS_py = neigh_Dc - targ_Dc 
        
        # compare IDL dLOS to dLOS_py 
        delta_dLOS = los_disp.dlos-dLOS_py
        print 'average difference in dLOS = ', np.mean(delta_dLOS)
        
        try: 
            combined_idl_dlos
        except NameError: 
            combined_idl_dlos = los_disp.dlos
            combined_py_dlos = dLOS_py
        else: 
            combined_idl_dlos = np.concatenate([combined_idl_dlos, los_disp.dlos]) 
            combined_py_dlos = np.concatenate([combined_py_dlos, dLOS_py])
        
        n_mocks = n_mocks+1
    
    #------------------------------------------------------------------------------------------
    # compare the dLOS histograms
    x_min = -1000.0
    x_max = 1000.0
    binsize = 0.2 
    n_bins = int((x_max-x_min)/binsize) 
    IDL_dlos_hist, mpc_binedges = np.histogram(combined_idl_dlos, bins=n_bins, range=[x_min, x_max])
    IDL_xlow = mpc_binedges[:-1]
    IDL_xhigh = mpc_binedges[1:] 
    IDL_xmid = np.array([0.5*(IDL_xlow[i]+IDL_xhigh[i]) for i in range(len(IDL_xlow))])
    
    PY_dlos_hist, mpc_binedges = np.histogram(combined_py_dlos, bins=n_bins, range=[x_min, x_max])
    PY_xlow = mpc_binedges[:-1]
    PY_xhigh = mpc_binedges[1:] 
    PY_xmid = np.array([0.5*(PY_xlow[i]+PY_xhigh[i]) for i in range(len(PY_xlow))])

    # plot dLOS ------------------------------------------------------------------------------
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    sub.plot(IDL_xmid, IDL_dlos_hist, lw=4, color=pretty_colors[-1], label=r"$d_{LOS}$ IDL") 
    sub.plot(PY_xmid, PY_dlos_hist, ls='--', lw=4, color=pretty_colors[0], label=r"$d_{LOS}$ COSMOLOPY") 

    sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
    sub.set_xlim([-20.0, 20.0])
    sub.set_ylim([0.0, 1.25*np.max(IDL_dlos_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_dlos_IDL_py_comparison.png', bbox_inches="tight")
    fig.clear() 

# ------------------------------------------------------------------------------------
# Derived quantities 
# ------------------------------------------------------------------------------------
def get_dlos_curvefit_sigma(fit='expon', binsize=0.1, **cat_corr): 
    ''' calculate sigma of the bestfit exponential/gaussian to dLOS histogram 
    dlos xrange is currently hardcoded
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    los_disp = dlos(**cat_corr) 
    los_disp.get_dlos_hist(binsize=binsize, normed=True) 
    dlos_xmid = los_disp.xmid 
    dlos_hist = los_disp.dlos_hist 
    
    if fit=='expon': 
        peak_xrange =  (dlos_xmid >= 0.0) & (dlos_xmid < 15.0)         # approximate peak xrange
        popt, pcov = curve_fit(los_disp.expon, dlos_xmid[peak_xrange], dlos_hist[peak_xrange])
    elif fit=='gauss':
        peak_xrange =  (dlos_xmid > -15.0) & (dlos_xmid < 15.0)         # approximate peak xrange
        popt, pcov = curve_fit(los_disp.gauss, dlos_xmid[peak_xrange], dlos_hist[peak_xrange])
    return popt[0]

def average_dlos_sigma(fit='expon', **cat_corr): 
    '''
    calculate average(sigma) of the bestfit exponential/gaussian of fiber-collided pair dLOS for mock catalogs
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    if catalog['name'].lower() == 'lasdamasgeo': 
        # For LasDamasGeo catalog 
        n_mock = 0
        sigma_sum = 0.0 
        for i_mock in range(1,41): 
            for letter in ['a', 'b', 'c', 'd']: 
                i_catalog = catalog.copy()
                i_catalog['n_mock'] = i_mock
                i_catalog['letter'] = letter
                i_cat_corr = {'catalog': i_catalog, 'correction': correction}  
                
                sigma_i = get_dlos_curvefit_sigma(fit=fit, **i_cat_corr) 
                
                n_mock = n_mock+1
                sigma_sum = sigma_sum + sigma_i

    return sigma_sum/np.float(n_mock) 

def average_dlos_fpeak(**cat_corr): 
    '''
    Calculate the average fpeak for certain catalog
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    f_1sigma = 0.0
    n_mock = 0
    if catalog['name'].lower() == 'lasdamasgeo': 
         for i_mock in range(1,41): 
             for letter in ['a', 'b', 'c', 'd']: 
                i_catalog = catalog.copy()
                i_catalog['n_mock'] = i_mock
                i_catalog['letter'] = letter
                i_cat_corr = {'catalog': i_catalog, 'correction': correction}  
                
                dlos_i = dlos(**i_cat_corr) 
                sigma_i = get_dlos_curvefit_sigma(**i_cat_corr) 
            
                # dlos within 1 sigma 
                dlos_i_1sigma = dlos_i.dlos[(dlos_i.dlos > -2.0*sigma_i) & (dlos_i.dlos < 2.0*sigma_i)]

                f_1sigma = f_1sigma+np.float(len(dlos_i_1sigma))/np.float(len(dlos_i.dlos))
                n_mock = n_mock+1
    else: 
        raise NameError('Not yet coded') 
    return f_1sigma/np.float(n_mock) 

def dlos_env(n=3, **cat_corr):
    ''' Using kth nearest neighbor distance compare dlos distrubtion  
    '''
    catalog = cat_corr['catalog']
    correction = {'name': 'upweight'} 
   
    # Read in dLOS data ----------------------------------------------------------------------
    if catalog['name'].lower() == 'qpm':    # QPM --------------------------------------------

        for i_mock in range(1, 11):         # loop through mocks 
            
            # read dLOS for each file 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
            print i_cat_corr

            los_disp_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
            # compute nth nearest neighbor distance for upweighted galaxy  
            NN_dist = genv.dlos_d_NN(n=n, **i_cat_corr ) 
            
            # combine dLOS and dNN from files 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i.dlos
                combined_dNN = NN_dist
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i.dlos])
                combined_dNN = np.concatenate([combined_dNN, NN_dist]) 

    else: 
        raise NotImplementedError('asdfasdfasdf') 

    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=(14,5))
    sub = fig.add_subplot(1,1,1)

    # loop through different dNN measurements 
    dNN_bins = [ (0, 5), (5, 10), (10, 50) ] 

    for i_bin, dNN_bin in enumerate(dNN_bins): 

        bin_index = ( combined_dNN >= dNN_bin[0] ) & ( combined_dNN < dNN_bin[1] ) 

        bin_dlos = combined_dlos[bin_index] 

        # Create histogram for combined dLOS values  (binsize is just guessed)
        x_min = -1000.0
        x_max = 1000.0
        binsize = 0.1 
        n_bins = int((x_max-x_min)/binsize) 
        
        # dLOS histograms 
        dlos_hist, mpc_binedges = np.histogram(bin_dlos, 
                bins=n_bins, range=[x_min, x_max], normed=True)

        mpc_low = mpc_binedges[:-1]
        mpc_high = mpc_binedges[1:] 
        mpc_mid = np.array([0.5*(mpc_low[i] + mpc_high[i]) for i in range(len(mpc_low))])
        
        sub.plot(mpc_mid, dlos_hist, 
                lw=4, color=pretty_colors[i_bin+1], 
                label='$'+str(dNN_bin[0])+' < d_'+str(n)+'NN < '+str(dNN_bin[1])+'$') 
        
    sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
    sub.set_xlim([-50.0, 50.0])
    sub.legend(loc='upper right') 

    fig.savefig(fc_util.get_fig_dir()+\
            'dlos_env_dependence_d'+str(n)+'NN.png', bbox_inches="tight")
    fig.clear()

#---- Fitting -----
def peak_expon(x, p): 
    '''
    Exponential function for the peak 
    '''
    p[0]*np.exp(-np.abs(x)/p[1])
    return p[0]*np.exp(-np.abs(x)/p[1])

def peak_gauss(x, p): 
    '''
    Exponential function for the peak 
    '''
    return p[0]*np.exp(-0.5*x**2/(p[1])**2)

def mpfit_peak_expon(p, fjac=None, x=None, y=None): 
    model = peak_expon(x, p) 
    status = 0 
    return([status, (y-model)]) 

def mpfit_peak_gauss(p, fjac=None, x=None, y=None): 
    model = peak_gauss(x, p) 
    status = 0 
    return([status, (y-model)]) 

def combined_dlos_fit(n_mocks, fit='gauss', sanitycheck=False, **cat_corr):
    '''
    Histogram using combined dLOS values from entire catalog, calculates appropriate binsize using Freedman-Diaconis rule, then fits the histogram
    Returns [sigma, fpeak]
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
   
    # Read in dLOS data ------------------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo ------------------------------------------------
        for i_mock in range(1,41): 
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

    elif catalog['name'].lower() == 'tilingmock':               # Tiling Mock ------------------------------------------------
        # import DLOS values for mock 
        los_disp = dlos(**cat_corr)

        combined_dlos = los_disp.dlos
        n_mocks = 1

    elif catalog['name'].lower() == 'qpm':                      # QPM ------------------------------------------------------------
        for i_mock in range(1,n_mocks+1): 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            los_disp_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
            # Combine dLOS values 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i.dlos
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i.dlos]) 

    elif catalog['name'].lower() == 'cmass': 
        los_disp = dlos(**cat_corr)
        combined_dlos = los_disp.dlos
        n_mocks = 1

    else: 
        raise NameError("not yet coded") 

    # --------------------------------------------------------------------------------------------------------------------------------
    # delete this later
    # write combined dlos 
    #np.savetxt(catalog['name'].lower()+str(n_mocks)+'mocks_combined_dlos.dat',            
    #        np.c_[combined_dlos], fmt=['%10.10f'], delimiter='\t')
    #------------------------------------------------------------------------------------------------------------------------------------ 
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
    #---------------------------------------------------------------------------------------------------------------------------------------
    # recompute histogram using freedman-diaconis binsize 
    n_bins = int((x_max-x_min)/binsize) 

    dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])

    # fitting ------------------------------------------------------------------------------------------------------------------------------------
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
    fpeak_xmin = -5.0*sigma
    fpeak_xmax = 5.0*sigma
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
    
    fpeak_dist = np.float(len(combined_dlos[(combined_dlos > -30.0) & (combined_dlos < 30.0)]))/np.float(len(combined_dlos))  # computed from the distribution

    print 'fpeak from distrib = ', fpeak_dist 

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
        #sub.set_xlim([-50.0, 50.0])
        sub.set_xlim([-5.0, 5.0])
        sub.set_ylim([0.0, 1.25*np.max(dlos_hist)])
        sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

        fig_dir = fc_util.get_fig_dir() 
        fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_dlos_peakfit_'+fit.lower()+'.png', bbox_inches="tight")
        fig.clear() 

    return [sigma, fpeak]

def combined_dlos_angsep_fit(n_mocks, fit='gauss', sanitycheck=False, **cat_corr):
    ''' Using combined dLOS values from entire catalog, we divide the dLOS values into bins of angular separation.  Then calculate appropriate binsize using Freedman-Diaconis rule, then fits the histogram
    Returns [sigma, fpeak]
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    angsep_cut = 0.01722/2.0
   
    # Read in dLOS data ------------------------------------------------------------------------------------------------------------
    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo ------------------------------------------------
        for i_mock in range(1,41): 
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

    elif catalog['name'].lower() == 'tilingmock':               # Tiling Mock ------------------------------------------------
        # import DLOS values for mock 
        los_disp = dlos(**cat_corr)

        combined_dlos = los_disp.dlos
        n_mocks = 1

    elif catalog['name'].lower() == 'qpm':                      # QPM ------------------------------------------------------------
        for i_mock in range(1,n_mocks+1): 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            los_disp_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
            
            # calculate the angular separate between the fiber-collided pairs 
            fc_angsep = ang_sep(los_disp_i.targ_ra, los_disp_i.targ_dec, los_disp_i.neigh_ra, los_disp_i.neigh_dec)

            # Combine dLOS values 
            try: 
                combined_dlos_close
                combined_dlos_far
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i.dlos
                combined_dlos_close = (los_disp_i.dlos)[fc_angsep < angsep_cut]
                combined_dlos_far = (los_disp_i.dlos)[fc_angsep >= angsep_cut]
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i.dlos])
                combined_dlos_close = np.concatenate([combined_dlos_close, (los_disp_i.dlos)[fc_angsep < angsep_cut]]) 
                combined_dlos_far = np.concatenate([combined_dlos_far, (los_disp_i.dlos)[fc_angsep >= angsep_cut]]) 

    elif catalog['name'].lower() == 'cmass': 
        los_disp = dlos(**cat_corr)
        combined_dlos = los_disp.dlos
        n_mocks = 1

    else: 
        raise NameError("not yet coded") 
    #------------------------------------------------------------------------------------------------------------------------------------ 
    # Determine appropriate binsize for dist (not too important)
    # Create histogram for combined dLOS values  (binsize is just guessed)
    x_min = -1000.0
    x_max = 1000.0
    binsize = 0.1 
    n_bins = int((x_max-x_min)/binsize) 
    
    # calculate the histograms 
    guess_dlos_hist_close, mpc_binedges = np.histogram(combined_dlos_close, bins=n_bins, range=[x_min, x_max])
    guess_dlos_hist_far, mpc_binedges = np.histogram(combined_dlos_far, bins=n_bins, range=[x_min, x_max])

    guess_xlow = mpc_binedges[:-1]
    guess_xhigh = mpc_binedges[1:] 
    guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
    
    # determine appropriate bin size using Freedman-Diaconis Rule 
    #peak_range = (guess_xmid >= 0.0) & (guess_xmid < 15.0) 
    #dlos_cumu = (guess_dlos_hist[peak_range]).cumsum() 
    #n_sample = dlos_cumu[-1]

    #iqr_index = fc_util.find_nearest(dlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    #iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
    #
    #binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
    binsize = 0.1
    print 'Freedman-Diaconis binsize = ', binsize
    #---------------------------------------------------------------------------------------------------------------------------------------
    # recompute histogram using freedman-diaconis binsize 
    n_bins = int((x_max-x_min)/binsize) 

    dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max], normed=True)
    dlos_hist_close, mpc_binedges = np.histogram(combined_dlos_close, bins=n_bins, range=[x_min, x_max], normed=True)
    dlos_hist_far, mpc_binedges = np.histogram(combined_dlos_far, bins=n_bins, range=[x_min, x_max], normed=True)
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])

    # fitting ------------------------------------------------------------------------------------------------------------------------------------
    #if fit.lower() == 'gauss': 
    #elif fit.lower() == 'expon': 
    #    peak_xrange = (xmid > -25.0) & (xmid < 25.0) 
    #else: 
    #    raise NameError('Error')
    #peak_xrange = (xmid > -25.0) & (xmid < 25.0)            # rough peak range 

    ## amplitude of dLOS distribution 
    #dlos_amp = np.mean(dlos_hist[(xmid >= -1.0) & (xmid < 1.0)])
    ##print 'Amplitude = ', dlos_amp, np.max(dlos_hist)
  
    ## MPFIT ------------------------------------------------------------------------------------
    #p0 = [dlos_amp, 5.0]            # initial guess

    #fa = {'x': xmid[peak_xrange], 'y': dlos_hist[peak_xrange]}
    #if fit.lower() == 'expon': 
    #    peak_pars = mpfit.mpfit(mpfit_peak_expon, p0, functkw=fa, nprint=0)
    #elif fit.lower() == 'gauss': 
    #    peak_pars = mpfit.mpfit(mpfit_peak_gauss, p0, functkw=fa, nprint=0)
    #else: 
    #    raise NameError("Fit not yet coded") 
    #
    #bestfit_amp = peak_pars.params[0]
    #sigma = peak_pars.params[1]
    #print fit.lower(), ' Best FIt Sigma = ', sigma
    #
    ##--------------------------------------------------------------------------------------------------------------------
    ## compute fpeak 
    #fpeak_xmin = -5.0*sigma
    #fpeak_xmax = 5.0*sigma
    #xrange = (xmid >= fpeak_xmin) & (xmid < fpeak_xmax) 
    ##fpeak = (np.sum(peak(xmid[xrange], popt[0])))/np.float(np.sum(dlos_hist)) 
    #if fit.lower() == 'expon': 
    #    fpeak = (np.sum(peak_expon(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    #elif fit.lower() == 'gauss': 
    #    fpeak = (np.sum(peak_gauss(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    #else: 
    #    raise NameError("Fit not yet coded") 
    #
    #fpeak_dist = np.float(len(combined_dlos[(combined_dlos > fpeak_xmin) & (combined_dlos < fpeak_xmax)]))/np.float(len(combined_dlos))  # computed from the distribution

    #print 'fpeak from fitting = ', fpeak
    #print 'fpeak from distrib = ', fpeak_dist 
    #
    #fpeak_dist = np.float(len(combined_dlos[(combined_dlos > -30.0) & (combined_dlos < 30.0)]))/np.float(len(combined_dlos))  # computed from the distribution

    #print 'fpeak from distrib = ', fpeak_dist 

    if sanitycheck == True:             # plot the distributions and fits for santicy check 
        prettyplot() 
        pretty_colors = prettycolors() 
        fig = plt.figure(1, figsize=[14,5]) 
        sub = fig.add_subplot(111) 

        ## fitting labels 
        #if fit.lower() == 'expon': 
        #    fit_label = "Exponential "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
        #else: 
        #    fit_label = "Gaussian "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
    
        #guess_scale = binsize/0.1
        #sub.plot(guess_xmid, guess_scale*guess_dlos_hist, lw=2, color=pretty_colors[-1], label=r"$d_{LOS}$ binsize $=0.1$") 

        sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[0], label=r"$d_{LOS}$ binsize $="+str(binsize)+"$ All")         # plot DLOS distribution
        sub.plot(xmid, dlos_hist_close, lw=4, color=pretty_colors[2], label=r"$d_{LOS}$ binsize $="+str(binsize)+"$ Ang Sep < 31")         # plot DLOS distribution
        sub.plot(xmid, dlos_hist_far, lw=4, color=pretty_colors[4], label=r"$d_{LOS}$ binsize $="+str(binsize)+"$ Ang Sep > 31")         # plot DLOS distribution

        # plot best fit 
        #if fit.lower() == 'expon': 
        #    sub.plot(xmid, peak_expon(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
        #elif fit.lower() == 'gauss': 
        #    sub.plot(xmid, peak_gauss(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
        #else: 
        #    raise NameError("Fit not yet coded") 

        # indicate the fpeak range (for sanity check purposes)
        #sub.vlines(fpeak_xmin, 0.0, 10.0*np.max(dlos_hist), color='k', lw=4)
        #sub.vlines(fpeak_xmax, 0.0, 10.0*np.max(dlos_hist), color='k', lw=4)

        #sub.text(-1.0*sigma, 0.25*np.max(dlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
        #sub.text(-1.0*sigma, 0.2*np.max(dlos_hist), r"$f_{peak,dist} = "+str(fpeak_dist)+"$") 
        sub.set_xlabel(r"$d_{LOS}$ (Mpc/h)", fontsize=20) 
        sub.set_xlim([-50.0, 50.0])
        #sub.set_xlim([-5.0, 5.0])
        sub.set_ylim([0.0, 1.25*np.max(dlos_hist_close)])
        sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

        fig_dir = fc_util.get_fig_dir() 
        fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_dlos_peakfit_'+fit.lower()+'_angsep_bins.png', bbox_inches="tight")
        fig.clear() 

    #return [sigma, fpeak]

def combined_dlos_dist(**cat_corr):
    # Combine dLOS values and output the normalized histogram for the peak of the distribution
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    n_mocks = 0 
    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo ------------------------------------------------
        for i_mock in range(1,41): 
            for letter in ['a', 'b', 'c', 'd']: 
                # individual catalog_correction dictonary 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock 
                i_catalog['letter'] = letter 
                i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
                
                los_d_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
                # Combine dLOS values 
                try: 
                    combined_dlos
                except NameError: 
                    combined_dlos = los_d_i.dlos
                else: 
                    combined_dlos = np.concatenate([combined_dlos, los_d_i.dlos]) 
                n_mocks = n_mocks+1

    elif catalog['name'].lower() == 'tilingmock':               # Tiling Mock ------------------------------------------------
        # import DLOS values for mock 
        los_d_i = dlos(**cat_corr)

        combined_dlos = los_d_i.dlos
        n_mocks = 1 

    elif catalog['name'].lower() == 'qpm':                      # QPM ------------------------------------------------------------
        for i_mock in range(1,101): 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            los_d_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
            # Combine dLOS values 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_d_i.dlos
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_d_i.dlos]) 
            n_mocks = n_mocks+1

    elif catalog['name'].lower() == 'cmass': 
        los_d_i = dlos(**cat_corr)
        combined_dlos = los_d_i.dlos

        n_mocks = 1
    else: 
        raise NameError("not yet coded") 
    
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
    peak_range = (guess_xmid >= 0.0) & (guess_xmid < 15.0)            # 5000 km/s is hardcoded for now 
    dlos_cumu = (guess_dlos_hist[peak_range]).cumsum() 
    n_sample = dlos_cumu[-1]

    iqr_index = fc_util.find_nearest(dlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
    
    binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
    print 'Freedman-Diaconis binsize = ', binsize
   
    # recompute histogram using appropriate binsize 
    # get normalized histogram ONLY near the peak
    new_x_min = -15.0
    new_x_max = 15.0
    narrow_range = (combined_dlos > new_x_min) & (combined_dlos < new_x_max)
    n_bins = int((new_x_max-new_x_min)/binsize) 

    print 'fpeak = ', np.float(len(combined_dlos[narrow_range]))/np.float(len(combined_dlos))

    dlos_hist, mpc_binedges = np.histogram(combined_dlos[narrow_range], bins=n_bins, range=[new_x_min, new_x_max], density=True)
    xlow = mpc_binedges[:-1]
    xhigh = mpc_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
    
    dlos_comb_peak_dist_file = ((los_d_i.file_name).rsplit('/',1))[0]+'/DLOS_norm_peak_dist_'+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined.dat'
    np.savetxt(dlos_comb_peak_dist_file,            
            np.c_[xmid, dlos_hist], fmt=['%10.5f', '%10.10f'], delimiter='\t')

def dlos_hist_binsize(**cat_corr): 
    '''
    Use Freedman-Diaconis rule to get the binsize for the peak distribution  
    '''
    # get line of sight dLOS histogram 
    LOSdisp = dlos(**cat_corr) 
    LOSdisp.get_dlos_hist(binsize=0.1) 
    LOSdisp_hist = LOSdisp.dlos_hist
    LOSdisp_xmid = LOSdisp.xmid
    
    # only consider what is "roughly" the peak (shouldn't matter too much because it tapers off fast) 
    peak_range = (LOSdisp_xmid >= 0.0) & (LOSdisp_xmid < 15.0) 
    LOSdisp_cumu = (LOSdisp_hist[peak_range]).cumsum() 

    # Need to determine the IQR 
    n_sample = LOSdisp_cumu[-1]
    iqr_index = fc_util.find_nearest(LOSdisp_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    iqr = 2.0*(LOSdisp_xmid[peak_range])[iqr_index]         #interquartile range 
    
    binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)
    return binsize 

def qpm_combined_dlos_zreal():
    ''' dLOS(z_real) distribution '''
    catalog = {'name': 'qpm'} 
    correction = {'name': 'true'} 
    
    n_mocks = 0 
    for i_mock in range(1, 11): 
        # individual catalog_correction dictonary 
        i_catalog = catalog.copy() 
        i_catalog['n_mock'] = i_mock 
        i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

        los_disp_rsd_i = dlos(readdata=False, **i_cat_corr)         # import DLOS values from each mock 
        los_disp_i_file = '.zreal.'.join((los_disp_rsd_i.file_name).rsplit('.', 1)) 

        if os.path.isfile(los_disp_i_file) == False: 
            build_dlos_zreal(**i_cat_corr) 
        
        los_disp_i = np.loadtxt(los_disp_i_file, unpack=True, usecols=[0]) 

        # Combine dLOS values 
        try: 
            combined_dlos
        except NameError: 
            combined_dlos = los_disp_i
        else: 
            combined_dlos = np.concatenate([combined_dlos, los_disp_i]) 
        n_mocks = n_mocks+1

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
   
    '''
    if fit.lower() == 'gauss': 
        peak_xrange = (xmid > -10.0) & (xmid < 10.0) 
    elif fit.lower() == 'expon': 
        peak_xrange = (xmid > -30.0) & (xmid < 30.0) 
    else: 
        raise NameError('Error')

    # UPDATED TO USE MPFIT
    # lambda functions for the fit  
    if fit.lower() == 'expon': 
        peak = lambda x, sig: dlos_amp*np.exp(-np.abs(x)/sig)
    elif fit.lower() == 'gauss': 
        peak = lambda x, sig: dlos_amp*np.exp(-0.5*x**2/sig**2)
    else: 
        raise NameError("Fit not yet coded") 
    
    # fit to dLOS histogram 
    popt, pcov = curve_fit(peak, xmid[peak_xrange], dlos_hist[peak_xrange]) #(updated to used MPFIT) 
    print 'SIGMA = ', popt[0]
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

    # compute fpeak 
    xrange = (xmid >= -15.0) & (xmid < 15.0) 
    #fpeak = (np.sum(peak(xmid[xrange], popt[0])))/np.float(np.sum(dlos_hist)) 
    if fit.lower() == 'expon': 
        fpeak = (np.sum(peak_expon(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    elif fit.lower() == 'gauss': 
        fpeak = (np.sum(peak_gauss(xmid[xrange], peak_pars.params)))/np.float(np.sum(dlos_hist)) 
    else: 
        raise NameError("Fit not yet coded") 
    
    # fitting labels 
    if fit.lower() == 'expon': 
        fit_label = "Exponential "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
    else: 
        fit_label = "Gaussian "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
    '''
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    guess_scale = binsize/0.1
    sub.plot(guess_xmid, guess_scale*guess_dlos_hist, lw=2, color=pretty_colors[-1], label=r"$d_{LOS}$ binsize $=0.1$") 
    sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[0], label=r"$d_{LOS}$ binsize $="+str(binsize)+"$") 

    '''
    if fit.lower() == 'expon': 
        sub.plot(xmid, peak_expon(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
    elif fit.lower() == 'gauss': 
        sub.plot(xmid, peak_gauss(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
    else: 
        raise NameError("Fit not yet coded") 
    '''

    #sub.text(-1.0*sigma, 0.25*np.max(dlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
    sub.set_xlabel(r"$d_{LOS}$", fontsize=20) 
    sub.set_xlim([-20.0, 20.0])
    sub.set_ylim([0.0, 1.25*np.max(dlos_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_dlos_zreal.png', bbox_inches="tight")
    fig.clear() 

def qpm_dlos_zbins(): 
    '''
    compare QPM dlos distribution over three redshift bins in order to determine if there are redshift dependence issues 
    everything hardcoded
    '''
    # set up figure 
    prettyplot() 
    pretty_colors = prettycolors() 

    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 
    
    # for low, mid, high zbin 
    for i_zbin, zbin in enumerate(['low', 'mid', 'high']): 

        # combine 10 mock dlos files
        n_mocks = 0 
        for i_file in range(1,11): 
            file_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/'        # directory
            dlos_file_name = ''.join([file_dir, 'dLOS_a0.6452_', str("%04d" % i_file), '.dr12d_cmass_ngc.vetoed.fibcoll.', zbin,'.dat']) 
            los_disp_i = np.loadtxt(dlos_file_name)

            # Combine dLOS values 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i]) 
            n_mocks = n_mocks+1

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
       
        # recompute histogram using appropriate binsize 
        n_bins = int((x_max-x_min)/binsize) 
        dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max], normed=True)
        print 'binsize ', binsize
        print n_bins
        xlow = mpc_binedges[:-1]
        xhigh = mpc_binedges[1:] 
        xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
        # amplitude of dLOS distribution 

        sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[i_zbin], label=r"$d_{LOS}$ "+zbin+" zbin") 

    #sub.text(-1.0*sigma, 0.25*np.max(dlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
    sub.set_xlabel(r"$d_{LOS}$", fontsize=20) 
    sub.set_xlim([-20.0, 20.0])
    sub.set_ylim([0.0, 1.25*np.max(dlos_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+'qpm_10mocks_combined_dlos_zbin_comparison.png', bbox_inches="tight")
    fig.clear() 

def ang_sep(ra1, dec1, ra2, dec2): 
    ''' Given a pair of ra and decs in DEGREES gives angular separation in DEGREES
    '''
    # convert to radians 
    ra1 = ra1*np.pi/180.
    dec1 = dec1*np.pi/180.
    ra2 = ra2*np.pi/180.
    dec2 = dec2*np.pi/180.

    x = np.cos(ra1)*np.cos(dec1)*np.cos(ra2)*np.cos(dec2) 
    y = np.sin(ra1)*np.cos(dec1)*np.sin(ra2)*np.cos(dec2) 
    z = np.sin(dec1)*np.sin(dec2)

    rad = np.arccos(x+y+z)
    
    sep = rad
    #sep = np.choose( rad<0.000004848 , ( np.sqrt( (np.cos(dec1)*(ra1-ra2))**2+(dec1-dec2)**2), rad))

    sep = sep*180./np.pi
    return sep

# ------------------------------------------------------------------------------------
# Plotting 
# ------------------------------------------------------------------------------------
def plot_dlos_peak(**cat_corr): 
    '''
    plot peak of dLOS distribution (testing appropriate binsize) 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    prettyplot() 
    pretty_colors = prettycolors()

    LOSdisp = dlos(**cat_corr) 
    LOSdisp.get_dlos_hist(binsize=0.25, normed=True) 
    sigma = get_dlos_curvefit_sigma(**cat_corr) 

    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    sub.plot(LOSdisp.xmid, LOSdisp.dlos_hist, color=pretty_colors[1]) 
    sub.plot(LOSdisp.xmid, LOSdisp.expon(LOSdisp.xmid, sigma), color=pretty_colors[5])
    sub.set_xlim([0.0, 5.0*sigma])  
    sub.set_ylim([0.0, 0.1])
    plt.show() 

def plot_dlos_tail(**cat_corr): 
    '''
    plot tail of dLOS distribution 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    prettyplot() 
    pretty_colors = prettycolors()

    LOSdisp = dlos(**cat_corr) 
    LOSdisp.get_dlos_hist(binsize=1.0, normed=True) 
    sigma = get_dlos_curvefit_sigma(**cat_corr) 

    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    sub.plot(LOSdisp.xmid, LOSdisp.dlos_hist, color=pretty_colors[1]) 
    sub.plot(LOSdisp.xmid, LOSdisp.expon(LOSdisp.xmid, sigma), color=pretty_colors[5])
    sub.set_xlim([0.0, 2.0*sigma])  
    sub.set_ylim([0.0, 0.2])
    plt.show() 

def plot_fcpaper_dlos(cat_corrs): 
    '''
    dLOS distribution of mocks for fc paper 
    Takes in multiple different files 
    '''
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=[14,5]) 
    sub = fig.add_subplot(111) 

    for icatcorr, cat_corr in enumerate(cat_corrs): 
        catalog = cat_corr['catalog']
        correction = cat_corr['correction']
        
        # number of mocks are hardcoded for the different mock catalogs 
        if catalog['name'].lower() == 'tilingmock': 
            mocks = ['']
        elif catalog['name'].lower() == 'lasdamasgeo': 
            mocks = [str(num)+letter for num in range(1,41) for letter in ['a', 'b', 'c', 'd']]
        elif catalog['name'].lower() == 'qpm': 
            mocks = range(1,101) 
        elif catalog['name'].lower() == 'cmass': 
            mocks = ['']
        else: 
            raise NameError('error')
    
        n_mocks = 0
        for i_mock in mocks: 
            i_catalog = catalog.copy()
            if catalog['name'].lower() != 'lasdamasgeo': 
                i_catalog['n_mock'] = i_mock 
            else: 
                i_catalog['n_mock'] = int(''.join(list(i_mock)[0:-1]))
                i_catalog['letter'] = list(i_mock)[-1]

            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
            
            los_disp_i = dlos(**i_cat_corr)         # import DLOS values from each mock 
            # Combine dLOS values 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i.dlos
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i.dlos]) 
            n_mocks = n_mocks+1

        # Create histogram for combined dLOS values  (binsize is just guessed)
        x_min = -1000.0
        x_max = 1000.0
        binsize = 0.2 
        n_bins = int((x_max-x_min)/binsize) 
    
        dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max], normed=True)
        xlow = mpc_binedges[:-1]
        xhigh = mpc_binedges[1:] 
        xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])

        lwid = 4
        if catalog['name'].lower() == 'lasdamasgeo': 
            cat_label = 'Las Damas'
            cat_color = pretty_colors[1]
        elif catalog['name'].lower() == 'qpm': 
            cat_label = 'QPM' 
            cat_color = pretty_colors[3]
        elif catalog['name'].lower() == 'tilingmock': 
            cat_label = 'Tiling Mock' 
            cat_color = pretty_colors[5]
        elif catalog['name'].lower() == 'cmass': 
            cat_label = 'BOSS CMASS'
            cat_color = pretty_colors[0]
            lwid = 2
        else: 
            raise NameError('asdf') 

        sub.plot(xmid, dlos_hist, lw=lwid, color=cat_color, label=cat_label) 

        del combined_dlos

    sub.set_xlabel(r"d$_{\rm{LOS}}$ (Mpc)", fontsize=24) 
    sub.set_xlim([-45.0, 45.0])
    #sub.set_xlim([0.1, 100.0])
    sub.set_ylim([0.0, 0.065])
    #sub.set_xscale('log') 
    #sub.set_yscale("log") 
    sub.legend(loc='upper left', scatterpoints=1, prop={'size':20}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+'fcpaper_dlos_dist.png', bbox_inches="tight")
    fig.clear() 

if __name__=="__main__": 
    cat_corr = {'catalog': {'name': 'qpm'}, 'correction': {'name': 'upweight'}} 
    dlos_env(n=1, **cat_corr)
    dlos_env(n=2, **cat_corr)
    dlos_env(n=3, **cat_corr)
    dlos_env(n=5, **cat_corr)

    correction = {'name': 'true'}
    cat_corrs = [{'catalog': {'name':catalog}, 'correction': correction} for catalog in ['cmass', 'lasdamasgeo', 'qpm', 'tilingmock']]

    cat_corrs = [{'catalog': {'name':catalog}, 'correction': correction} for catalog in ['qpm']]
    #for cat_corr in cat_corrs: 
    #    combined_dlos_dist(**cat_corr)

    #plot_fcpaper_dlos(cat_corrs)

    #print 'LASDAMASGEO------------------------------------------------------'
    #cat_corr = {'catalog': {'name':'lasdamasgeo'}, 'correction': {'name': 'upweight'}}
    #print 'Expon ', combined_dlos_fit(40, fit='expon', sanitycheck=True,  **cat_corr) 
    #print 'Gauss ', combined_dlos_fit(40, fit='gauss', sanitycheck=True,  **cat_corr) 
    print 'QPM------------------------------------------------------'
    cat_corr = {'catalog': {'name':'qpm'}, 'correction': {'name': 'upweight'}}
    #print 'Expon ', combined_dlos_fit(10, fit='expon', sanitycheck=True,  **cat_corr) 
    #print 'Expon ', combined_dlos_angsep_fit(10, fit='expon', sanitycheck=True,  **cat_corr) 
    #print 'Gauss ', combined_dlos_fit(100, fit='gauss', sanitycheck=True,  **cat_corr) 
    #print 'Tiling Mock------------------------------------------------------'
    #cat_corr = {'catalog': {'name':'tilingmock'}, 'correction': {'name': 'upweight'}}
    #print 'Expon ', combined_dlos_fit(1, fit='expon', sanitycheck=True,  **cat_corr) 
    #print 'Gauss ', combined_dlos_fit(1, fit='gauss', sanitycheck=True,  **cat_corr) 
    #cat_corr = {'catalog': {'name':'cmass'}, 'correction': {'name': 'upweight'}}
    #print 'CMASS-----------------------------------------------------------'
    #print 'Expon ', combined_dlos_fit(1, fit='expon', sanitycheck=True,  **cat_corr) 
    #print 'Gauss ', combined_dlos_fit(1, fit='gauss', sanitycheck=True,  **cat_corr) 
