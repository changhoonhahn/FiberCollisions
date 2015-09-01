import numpy as np
from scipy.optimize import curve_fit
import sys
import os.path
import subprocess
import cosmolopy as cosmos
import fibcol_data as fc_data
import fibcol_utility as fc_util
import mpfit as mpfit

# ------------------------------------------------------------------------------------
# Set up LOS displacement class
class vlos:
    def __init__(self, readdata=True, vlos_method='cDz', **cat_corr): 
        '''
        Given catalog_correction info, read vLOS values. If vLOS values do not exist, make them
        '''
        catalog = cat_corr['catalog']
        correction = cat_corr['correction']

        # LasDamas Geo ------------------------------------------------------------------------------------------
        if catalog['name'].lower() == 'lasdamasgeo': 

            file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
            #File name 
            file_name = file_dir+'vLOS_sdssmock_gamma_lrgFull_zm_oriana'+str(catalog['n_mock']+100)[1:3]+catalog['letter']+'_no.rdcz.dat'
            self.file_name = file_name 
            
            if readdata == True: 
                if os.path.isfile(file_name) == False: 
                    # if file does not exist then
                    build_vlos(**cat_corr) 

                self.vlos = np.loadtxt(file_name) 
        
        # Tiling Mock ------------------------------------------------------------------------------------------
        elif catalog['name'].lower() == 'tilingmock': 

            file_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'        # directory
            # File name 
            file_name = file_dir+'vLOS_cmass-boss5003sector-icoll012.dat'
            self.file_name = file_name 
            
            if readdata == True: 
                if os.path.isfile(file_name) == False: 
                    # if file does not exist then
                    build_vlos(**cat_corr) 

                self.vlos = np.loadtxt(file_name) 
        
        # QPM ----------------------------------------------------------------------------------------------------
        elif catalog['name'].lower() == 'qpm': 
            file_dir = '/mount/riachuelo1/hahn/data/QPM/dr12d/'        # directory
            # File name 
            file_name = ''.join([file_dir, 'vLOS_a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.vetoed.fibcoll.dat']) 
            self.file_name = file_name 
            
            if readdata == True: 
                if os.path.isfile(file_name) == False: 
                    print 'Constructing ', file_name 
                    # if file does not exist then
                    build_vlos(**cat_corr) 
                
                # to account for the different methods of calculating vLOS 
                if vlos_method == 'cDz': 
                    self.vlos = np.loadtxt(file_name, unpack=True, usecols=[0])
                elif vlos_method == 'z2v': 
                    self.vlos = np.loadtxt(file_name, unpack=True, usecols=[1])
                elif vlos_method == 'vpec': 
                    self.vlos = np.loadtxt(file_name, unpack=True, usecols=[2])
                else: 
                    raise NameError('wtf is going on') 
        
        # CMASS --------------------------------------------------------------------------------------------------------
        elif catalog['name'].lower() == 'cmass': 
            file_dir = '/mount/riachuelo1/hahn/data/'        # directory
            # File name 
            file_name = ''.join([file_dir, 'vLOS_cmass-dr12v4-N-Reid-weights-zlim.dat']) 
            self.file_name = file_name 
            
            if readdata == True: 
                if os.path.isfile(file_name) == False: 
                    print 'Constructing ', file_name 
                    # if file does not exist then
                    build_vlos(**cat_corr) 

                self.vlos = np.loadtxt(file_name) 

        else: 
            raise NameError('not coded yet') 

    def get_vlos_hist(self, binsize=1000.0, normed=False): 
        ''' 
        Compute the histogram for vlos  class
        '''
        x_min = -80000.0
        x_max = 80000.0 
        n_bins = int((x_max-x_min)/binsize) 
        vlos_hist, mpc_binedges = np.histogram(self.vlos, bins=n_bins, range=[x_min, x_max], normed=normed)
        # store histogram info 
        self.vlos_hist = vlos_hist 
        self.xlow = mpc_binedges[:-1]
        self.xhigh = mpc_binedges[1:] 
        self.xmid = np.array([0.5*(self.xlow[i]+self.xhigh[i]) for i in range(len(self.xlow))])

    def expon(self, x, sig):
        return np.max(self.vlos_hist)*np.exp(-np.abs(x)/sig)
    def gauss(self, x, sig):
        return np.max(self.vlos_hist)*np.exp(-0.5*x**2/sig**2)

def build_vlos(**cat_corr): 
    '''
    build DLOS for given catalog_correction ID 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    mock = fc_data.galaxy_data('data', readdata=True, **cat_corr)
    mock_file = mock.file_name 
    los_disp = vlos(readdata=False, **cat_corr)  
    vlos_file = los_disp.file_name

    build_vlos_idl_cmd = ''. join(['idl -e ', '"', 'build_fibcoll_vlos, ', "'", catalog['name'], "','", mock_file, "','", vlos_file, "'", '"'])
    print build_vlos_idl_cmd 
    os.system(build_vlos_idl_cmd) 

def build_vlos_zreal(**cat_corr): 
    '''
    build DLOS for given catalog_correction ID 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    if catalog['name'].lower() != 'qpm': raise NameError('asdfasdfasdf') 

    mock_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
        'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz']) 

    los_disp = vlos(readdata=False, **cat_corr)  
    vlos_file = '.zreal.'.join((los_disp.file_name).rsplit('.', 1)) 

    build_vlos_idl_cmd = ''. join(['idl -e ', '"', 'build_fibcoll_vlos_zreal, ', "'", catalog['name'], "','", mock_file, "','", vlos_file, "'", '"'])
    print build_vlos_idl_cmd 
    os.system(build_vlos_idl_cmd) 

def build_vlos_zreal_comp(**cat_corr): 
    '''
    build DLOS for given catalog_correction ID 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    if catalog['name'].lower() != 'qpm': raise NameError('asdfasdfasdf') 

    mock_file = ''.join(['/mount/riachuelo2/rs123/BOSS/QPM/cmass/mocks/dr12d/ngc/data/', 
        'a0.6452_', str("%04d" % catalog['n_mock']), '.dr12d_cmass_ngc.rdz']) 

    los_disp = vlos(readdata=False, **cat_corr)  
    vlos_file = '.zreal.comp.'.join((los_disp.file_name).rsplit('.', 1)) 

    build_vlos_idl_cmd = ''. join(['idl -e ', '"', 'build_fibcoll_vlos_zreal_comp, ', "'", catalog['name'], "','", mock_file, "','", vlos_file, "'", '"'])
    print build_vlos_idl_cmd 
    os.system(build_vlos_idl_cmd) 

# ------------------------------------------------------------------------------------
# Derived quantities 
def get_vlos_curvefit_sigma(fit='expon', binsize=100.0, **cat_corr): 
    '''
    calculate sigma of the bestfit exponential/gaussian to vLOS histogram 
    vLOS xrange is currently hardcoded
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    los_disp = vlos(**cat_corr) 
    los_disp.get_vlos_hist(binsize=binsize, normed=True) 
    vlos_xmid = los_disp.xmid 
    vlos_hist = los_disp.vlos_hist 
    
    if fit=='expon': 
        peak_xrange =  (vlos_xmid >= 0.0) & (vlos_xmid < 30000.0)         # approximate peak xrange
        popt, pcov = curve_fit(los_disp.expon, vlos_xmid[peak_xrange], vlos_hist[peak_xrange])
    elif fit=='gauss':
        peak_xrange =  (vlos_xmid > -30000.0) & (vlos_xmid < 30000.0)         # approximate peak xrange
        popt, pcov = curve_fit(los_disp.gauss, vlos_xmid[peak_xrange], vlos_hist[peak_xrange])
    return popt[0]

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

def peak_loren(x, p): 
    '''
    Exponential function for the peak 
    '''
    return p[0]/(1.0+(x/p[1])**2)

def mpfit_peak_expon(p, fjac=None, x=None, y=None): 
    model = peak_expon(x, p) 
    status = 0 
    return([status, (y-model)]) 

def mpfit_peak_gauss(p, fjac=None, x=None, y=None): 
    model = peak_gauss(x, p) 
    status = 0 
    return([status, (y-model)]) 

def mpfit_peak_loren(p, fjac=None, x=None, y=None): 
    model = peak_loren(x, p) 
    status = 0 
    return([status, (y-model)]) 

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

def combined_vlos_fit(fit='gauss', sanitycheck=False, **cat_corr):
    '''
    Histogram using combined vLOS values from entire catalog, calculates appropriate binsize using Freedman-Diaconis rule, then fits the histogram
    Returns [sigma, fpeak]
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    peak_min = -3000.0
    peak_max = 3000.0

    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1, figsize=[14,5]) 
    sub = fig.add_subplot(111) 
    # for i_vlos, vlos_method in enumerate(['cDz', 'vpec']):          # currently hardcoded to do calculations for both cDz and vpec

    n_mocks = 0 
    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo ------------------------------------------------
        for i_mock in range(1,41): 
            for letter in ['a', 'b', 'c', 'd']: 
                # individual catalog_correction dictonary 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock 
                i_catalog['letter'] = letter 
                i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
                
                los_v_i = vlos(**i_cat_corr)         # import DLOS values from each mock 
                # Combine dLOS values 
                try: 
                    combined_vlos
                except NameError: 
                    combined_vlos = los_v_i.vlos
                else: 
                    combined_vlos = np.concatenate([combined_vlos, los_v_i.vlos]) 
                n_mocks = n_mocks+1

    elif catalog['name'].lower() == 'tilingmock':               # Tiling Mock ------------------------------------------------
        # import DLOS values for mock 
        los_v = vlos(**cat_corr)

        combined_vlos = los_v.vlos
        n_mocks = n_mocks+1

    elif catalog['name'].lower() == 'qpm':                      # QPM ------------------------------------------------------------
        for i_mock in range(1,101): 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            #los_v_i = vlos(vlos_method=vlos_method, **i_cat_corr)         # import DLOS values from each mock 
            los_v_i = vlos(vlos_method='cDz', **i_cat_corr)         # import vLOS values from each mock (cDz)

            # Combine vLOS values 
            try: 
                combined_vlos
            except NameError: 
                combined_vlos = los_v_i.vlos
            else: 
                combined_vlos = np.concatenate([combined_vlos, los_v_i.vlos]) 
            n_mocks = n_mocks+1

    elif catalog['name'].lower() == 'cmass': 
        los_v = vlos(**cat_corr)
        combined_vlos = los_v.vlos
        n_mocks = n_mocks+1

    else: 
        raise NameError("not yet coded") 
    
    # Create histogram for combined dLOS values  (binsize is just guessed)
    x_min = -80000.0
    x_max = 80000.0
    binsize = 10.0 
    n_bins = int((x_max-x_min)/binsize) 
    
    guess_vlos_hist, v_binedges = np.histogram(combined_vlos, bins=n_bins, range=[x_min, x_max])
    guess_xlow = v_binedges[:-1]
    guess_xhigh = v_binedges[1:] 
    guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
    
    # determine appropriate bin size using Freedman-Diaconis Rule 
    peak_range = (guess_xmid >= 0.0) & (guess_xmid < 5000.0)            # 5000 km/s is hardcoded for now 
    vlos_cumu = (guess_vlos_hist[peak_range]).cumsum() 
    n_sample = vlos_cumu[-1]

    iqr_index = fc_util.find_nearest(vlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
    
    binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
    if binsize < 30.0: 
        binsize = 30.0
    print 'Freedman-Diaconis binsize = ', binsize
   
    # recompute histogram using appropriate binsize 
    new_x_min = -5000.0
    new_x_max = 5000.0
    narrow_range = (combined_vlos > new_x_min) & (combined_vlos < new_x_max) 
    n_bins = int((new_x_max-new_x_min)/binsize) 
    vlos_hist, v_binedges = np.histogram(combined_vlos[narrow_range], bins=n_bins, range=[new_x_min, new_x_max], normed=True)
    print 'binsize', binsize
    print 'nbin', int((x_max-x_min)/binsize) 
    #print v_binedges
    xlow = v_binedges[:-1]
    xhigh = v_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])

    peak_xrange = (xmid > peak_min) & (xmid < peak_max) 

    # amplitude of dLOS distribution 
    vlos_amp = np.mean(vlos_hist[(xmid >= -100.0) & (xmid < 100.0)])
    print vlos_amp, np.max(vlos_hist)
    
    p0 = [vlos_amp, 500.0]            # initial guess

    fa = {'x': xmid[peak_xrange], 'y': vlos_hist[peak_xrange]}
    if fit.lower() == 'expon': 
        peak_pars = mpfit.mpfit(mpfit_peak_expon, p0, functkw=fa, nprint=0)
    elif fit.lower() == 'gauss': 
        peak_pars = mpfit.mpfit(mpfit_peak_gauss, p0, functkw=fa, nprint=0)
    elif fit.lower() == 'loren': 
        peak_pars = mpfit.mpfit(mpfit_peak_loren, p0, functkw=fa, nprint=0)
    else: 
        raise NameError("Fit not yet coded") 
    
    bestfit_amp = peak_pars.params[0]
    sigma = peak_pars.params[1]
    print fit.lower()
    print 'bestfit sigma = ', sigma

    # compute fpeak 
    xrange = (xmid >= -5.0*sigma) & (xmid < 5.0*sigma) 
    if fit.lower() == 'expon': 
        peak_sum = np.sum(peak_expon(xmid[xrange], peak_pars.params)) 
    elif fit.lower() == 'gauss': 
        peak_sum = np.sum(peak_gauss(xmid[xrange], peak_pars.params)) 
    elif fit.lower() == 'loren': 
        peak_sum = np.sum(peak_loren(xmid[xrange], peak_pars.params)) 
    else: 
        raise NameError("Fit not yet coded") 

    fpeak = peak_sum/np.sum(vlos_hist)

    print '+/-5 sigma x-range ', -5.0*sigma, '< v < ', 5.0*sigma
    print 'fpeak, calculated from fit = ', fpeak
    fpeak_dist = np.float(len(combined_vlos[(combined_vlos > -5.0*sigma) & (combined_vlos < 5.0*sigma)]))/np.float(len(combined_vlos))
    print 'fpeak, calculated from dist = ', fpeak_dist 

    if np.abs(fpeak_dist/fpeak - 1.0) > 0.5: 
        print 'fpeak values do not agree very well'

    # compute the chi-squared value for peak range
    if fit.lower() == 'expon': 
        E_peak = peak_expon(xmid, peak_pars.params)
    elif fit.lower() == 'gauss': 
        E_peak = peak_gauss(xmid, peak_pars.params)
    elif fit.lower() == 'loren': 
        E_peak = peak_gauss(xmid, peak_pars.params)
    else: 
        raise NameError("Fit not yet coded") 

    #print vlos_hist[peak_xrange]
    #print E_peak[peak_xrange]
    #print (vlos_hist[peak_xrange]-E_peak[peak_xrange])**2/E_peak[peak_xrange]
    chi2 = 1.0/np.float(len(xmid[peak_xrange])-3)*(np.sum((vlos_hist[peak_xrange]-E_peak[peak_xrange])**2 / E_peak[peak_xrange]))
    #print 'Reduced Chi-Squared = ', chi2 

    if sanitycheck == True: 
    
        sub.plot(xmid, vlos_hist, lw=4, color=pretty_colors[0], label=r"$v_{LOS}$ binsize $="+str(binsize)+"$") 
        #sub.plot(xmid, vlos_hist, lw=4, color=pretty_colors[i_vlos], label=r"$v_{LOS}$ "+vlos_method) 
        
        # plot best fit function  
        if fit.lower() == 'expon': 
            fit_label = "Exponential "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
            sub.plot(xmid, peak_expon(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)

        elif fit.lower() == 'gauss': 
            fit_label = "Gaussian "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
            sub.plot(xmid, peak_gauss(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)

        elif fit.lower() == 'loren': 
            fit_label = "Lorentzian "+r"$(\gamma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
            sub.plot(xmid, peak_loren(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)

        else: 
            raise NameError("Fit not yet coded") 
    
        #sub.text(-1750.0, 0.5*np.max(vlos_hist), 'peak: '+str(peak_min)+' - '+str(peak_max))
        #sub.text(-1750.0, 0.45*np.max(vlos_hist), 'reduced chi2: '+str(chi2)) 
        if fit.lower() != 'loren': 
            sub.text(-1.0*sigma, 0.25*np.max(vlos_hist), r"$f_{peak} = "+str(fpeak_dist)+"$") 
        else: 
            sub.text(-500.0, 0.25*np.max(vlos_hist), r"$f_{peak} = "+str(fpeak_dist)+"$") 
        sub.vlines(-5.0*sigma, 0, 1000000, colors='k', lw='4')
        sub.vlines(5.0*sigma, 0, 1000000, colors='k', lw='4')

    sub.set_xlabel(r"$v_{LOS}$", fontsize=20) 
    sub.set_xlim([-7500.0, 7500.0])
    sub.set_ylim([0.0, 1.25*np.max(vlos_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_dir = fc_util.get_fig_dir() 
    #fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_vlos_comp_peakfit_'+fit.lower()+'.png', bbox_inches="tight")
    #fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_vlos_'+vlos_method+'_peakfit_'+fit.lower()+'.png', bbox_inches="tight")
    fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_vlos_peakfit_'+fit.lower()+'.png', bbox_inches="tight")
    fig.clear() 

    return [sigma, fpeak]

def combined_vlos_dist(**cat_corr):
    '''
    Using combined vLOS values, outputs the normalized histogram for the peak of the distribution to be used for vlos sampling 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    peak_min = -2000.0
    peak_max = 2000.0

    n_mocks = 0 
    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo ------------------------------------------------
        for i_mock in range(1,41): 
            for letter in ['a', 'b', 'c', 'd']: 
                # individual catalog_correction dictonary 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock 
                i_catalog['letter'] = letter 
                i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
                
                los_v_i = vlos(**i_cat_corr)         # import DLOS values from each mock 
                # Combine dLOS values 
                try: 
                    combined_vlos
                except NameError: 
                    combined_vlos = los_v_i.vlos
                else: 
                    combined_vlos = np.concatenate([combined_vlos, los_v_i.vlos]) 
                n_mocks = n_mocks+1

    elif catalog['name'].lower() == 'tilingmock':               # Tiling Mock ------------------------------------------------
        # import DLOS values for mock 
        los_v = vlos(**cat_corr)

        combined_vlos = los_v.vlos
        n_mocks = n_mocks+1

    elif catalog['name'].lower() == 'qpm':                      # QPM ------------------------------------------------------------
        for i_mock in range(1,101): 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            los_v_i = vlos(**i_cat_corr)         # import DLOS values from each mock 
            # Combine dLOS values 
            try: 
                combined_vlos
            except NameError: 
                combined_vlos = los_v_i.vlos
            else: 
                combined_vlos = np.concatenate([combined_vlos, los_v_i.vlos]) 
            n_mocks = n_mocks+1

    elif catalog['name'].lower() == 'cmass': 
        los_v_i = vlos(**cat_corr)
        combined_vlos = los_v_i.vlos
        n_mocks = n_mocks+1

    else: 
        raise NameError("not yet coded") 
    
    # Create histogram for combined dLOS values  (binsize is just guessed)
    x_min = -80000.0
    x_max = 80000.0
    binsize = 10.0 
    n_bins = int((x_max-x_min)/binsize) 
    
    guess_vlos_hist, v_binedges = np.histogram(combined_vlos, bins=n_bins, range=[x_min, x_max])
    guess_xlow = v_binedges[:-1]
    guess_xhigh = v_binedges[1:] 
    guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
    
    # determine appropriate bin size using Freedman-Diaconis Rule 
    peak_range = (guess_xmid >= 0.0) & (guess_xmid < 5000.0)            # 5000 km/s is hardcoded for now 
    vlos_cumu = (guess_vlos_hist[peak_range]).cumsum() 
    n_sample = vlos_cumu[-1]

    iqr_index = fc_util.find_nearest(vlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
    
    binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
    if binsize < 30.0: 
        binsize = 30.0
    print 'Freedman-Diaconis binsize = ', binsize
   
    # recompute histogram using appropriate binsize 
    new_x_min = -5000.0
    new_x_max = 5000.0
    narrow_range = (combined_vlos > new_x_min) & (combined_vlos < new_x_max) 
    n_bins = int((new_x_max-new_x_min)/binsize) 

    vlos_hist, v_binedges = np.histogram(combined_vlos[narrow_range], bins=n_bins, range=[new_x_min, new_x_max], density=True)
    vlos_hist = vlos_hist*binsize       # so that the integral adds up to 1.0
    
    xlow = v_binedges[:-1]
    xhigh = v_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
    
    vlos_comb_peak_dist_file = ((los_v_i.file_name).rsplit('/',1))[0]+'/vLOS_norm_peak_dist_'+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined.dat'
    np.savetxt(vlos_comb_peak_dist_file,            
            np.c_[xmid, vlos_hist], fmt=['%10.5f', '%10.10f'], delimiter='\t')

def qpm_combined_vlos_zreal():
    ''' '''
    catalog = {'name': 'qpm'} 
    correction = {'name': 'true'} 
    
    n_mocks = 0 
    for i_mock in range(1, 11): 
        # individual catalog_correction dictonary 
        i_catalog = catalog.copy() 
        i_catalog['n_mock'] = i_mock 
        i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

        los_disp_rsd_i = vlos(readdata=False, **i_cat_corr)         # import DLOS values from each mock 
        los_disp_i_file = '.zreal.'.join((los_disp_rsd_i.file_name).rsplit('.', 1)) 

        if os.path.isfile(los_disp_i_file) == False: 
            build_vlos_zreal(**i_cat_corr) 
        
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
    x_min = -80000.0
    x_max = 80000.0
    binsize = 10.0 
    n_bins = int((x_max-x_min)/binsize) 
    
    guess_dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
    guess_xlow = mpc_binedges[:-1]
    guess_xhigh = mpc_binedges[1:] 
    guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
    
    # determine appropriate bin size using Freedman-Diaconis Rule 
    peak_range = (guess_xmid >= 0.0) & (guess_xmid < 5000.0) 
    dlos_cumu = (guess_dlos_hist[peak_range]).cumsum() 
    n_sample = dlos_cumu[-1]

    iqr_index = fc_util.find_nearest(dlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
    iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
    
    binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
    if binsize < 30.0: 
        binsize = 30.0
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
   
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    guess_scale = binsize/10.0
    sub.plot(guess_xmid, guess_scale*guess_dlos_hist, lw=2, color=pretty_colors[-1], label=r"$v_{LOS}$ binsize $=0.1$") 
    sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[0], label=r"$v_{LOS}$ binsize $="+str(binsize)+"$") 

    #sub.text(-1.0*sigma, 0.25*np.max(dlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
    sub.set_xlabel(r"$v_{LOS}$", fontsize=20) 
    sub.set_xlim([-3000.0, 3000.0])
    sub.set_ylim([0.0, 1.25*np.max(dlos_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_vlos_zreal.png', bbox_inches="tight")
    fig.clear() 

def qpm_combined_vlos_zreal_comp():
    ''' 
    QPM comparison vLOS(zreal) using three different methods to compute vLOS 
    
    '''
    catalog = {'name': 'qpm'} 
    correction = {'name': 'true'} 
        
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 
   
    for i_vlos, vlos_method in enumerate(['c Dz', 'z2v', 'vpec']): 
        n_mocks = 0 
        for i_mock in range(1, 11): 
            # individual catalog_correction dictonary 
            i_catalog = catalog.copy() 
            i_catalog['n_mock'] = i_mock 
            i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

            los_disp_rsd_i = vlos(readdata=False, **i_cat_corr)         # import DLOS values from each mock 
            los_disp_i_file = '.zreal.comp.'.join((los_disp_rsd_i.file_name).rsplit('.', 1)) 

            if os.path.isfile(los_disp_i_file) == False: 
                build_vlos_zreal_comp(**i_cat_corr) 
            
            los_disp_i = np.loadtxt(los_disp_i_file, unpack=True, usecols=[i_vlos]) 

            # Combine dLOS values 
            try: 
                combined_dlos
            except NameError: 
                combined_dlos = los_disp_i
            else: 
                combined_dlos = np.concatenate([combined_dlos, los_disp_i]) 
            n_mocks = n_mocks+1

        # Create histogram for combined dLOS values  (binsize is just guessed)
        x_min = -80000.0
        x_max = 80000.0
        binsize = 10.0 
        n_bins = int((x_max-x_min)/binsize) 
        
        guess_dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max])
        guess_xlow = mpc_binedges[:-1]
        guess_xhigh = mpc_binedges[1:] 
        guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
        
        # determine appropriate bin size using Freedman-Diaconis Rule 
        peak_range = (guess_xmid >= 0.0) & (guess_xmid < 5000.0) 
        dlos_cumu = (guess_dlos_hist[peak_range]).cumsum() 
        n_sample = dlos_cumu[-1]

        iqr_index = fc_util.find_nearest(dlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
        iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
        
        binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
        if binsize < 30.0: 
            binsize = 30.0
        print 'Freedman-Diaconis binsize = ', binsize
       
        # recompute histogram using appropriate binsize 
        n_bins = int((x_max-x_min)/binsize) 
        dlos_hist, mpc_binedges = np.histogram(combined_dlos, bins=n_bins, range=[x_min, x_max], normed=True)
        print 'binsize ', binsize
        print n_bins
        xlow = mpc_binedges[:-1]
        xhigh = mpc_binedges[1:] 
        xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
        
        if vlos_method != 'z2v': 
            # amplitude of dLOS distribution 
            sub.plot(xmid, dlos_hist, lw=4, color=pretty_colors[i_vlos], label=r"$v_{LOS}$ "+vlos_method) 

    #sub.text(-1.0*sigma, 0.25*np.max(dlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
    sub.set_xlabel(r"$v_{LOS}$", fontsize=20) 
    sub.set_xlim([-3000.0, 3000.0])
    sub.set_ylim([0.0, 1.25*np.max(dlos_hist)])
    sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_vlos_zreal_comp.png', bbox_inches="tight")
    fig.clear() 

def combined_vlos_fit_peakrange_test(fit='gauss', sanitycheck=False, **cat_corr):
    '''
    Histogram using combined vLOS values from entire catalog, calculates appropriate binsize using Freedman-Diaconis rule, then fits the histogram
    Returns [sigma, fpeak]
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']

    for i_peak in range(1, 4): 
        fig = plt.figure(1, figsize=[15,15]) 
        peak_min = -1.0*(1.0+0.5*np.float(i_peak-1))*1000.0
        peak_max = 1.0*(1.0+0.5*np.float(i_peak-1))*1000.0
        print peak_min, peak_max

        n_mocks = 0 
        if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo ------------------------------------------------
            for i_mock in range(1,41): 
                for letter in ['a', 'b', 'c', 'd']: 
                    # individual catalog_correction dictonary 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock 
                    i_catalog['letter'] = letter 
                    i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
                    
                    los_v_i = vlos(**i_cat_corr)         # import DLOS values from each mock 
                    # Combine dLOS values 
                    try: 
                        combined_vlos
                    except NameError: 
                        combined_vlos = los_v_i.vlos
                    else: 
                        combined_vlos = np.concatenate([combined_vlos, los_v_i.vlos]) 
                    n_mocks = n_mocks+1

        elif catalog['name'].lower() == 'tilingmock':               # Tiling Mock ------------------------------------------------
            # import DLOS values for mock 
            los_v = vlos(**cat_corr)

            combined_vlos = los_v.vlos
            n_mocks = n_mocks+1

        elif catalog['name'].lower() == 'qpm':                      # QPM ------------------------------------------------------------
            for i_mock in range(1,101): 
                # individual catalog_correction dictonary 
                i_catalog = catalog.copy() 
                i_catalog['n_mock'] = i_mock 
                i_cat_corr = {'catalog':i_catalog, 'correction': correction} 

                los_v_i = vlos(**i_cat_corr)         # import DLOS values from each mock 
                # Combine dLOS values 
                try: 
                    combined_vlos
                except NameError: 
                    combined_vlos = los_v_i.vlos
                else: 
                    combined_vlos = np.concatenate([combined_vlos, los_v_i.vlos]) 
                n_mocks = n_mocks+1

        elif catalog['name'].lower() == 'cmass': 
            los_v = vlos(**cat_corr)
            combined_vlos = los_v.vlos
            n_mocks = n_mocks+1

        else: 
            raise NameError("not yet coded") 
        
        # Create histogram for combined dLOS values  (binsize is just guessed)
        x_min = -50000.0
        x_max = 50000.0
        binsize = 10.0 
        n_bins = int((x_max-x_min)/binsize) 
        
        guess_vlos_hist, v_binedges = np.histogram(combined_vlos, bins=n_bins, range=[x_min, x_max])
        guess_xlow = v_binedges[:-1]
        guess_xhigh = v_binedges[1:] 
        guess_xmid = np.array([0.5*(guess_xlow[i]+guess_xhigh[i]) for i in range(len(guess_xlow))])
        
        # determine appropriate bin size using Freedman-Diaconis Rule 
        peak_range = (guess_xmid >= 0.0) & (guess_xmid < peak_max) 
        vlos_cumu = (guess_vlos_hist[peak_range]).cumsum() 
        n_sample = vlos_cumu[-1]

        iqr_index = fc_util.find_nearest(vlos_cumu, np.int(np.floor(n_sample/2.0)), index=True)
        iqr = 2.0*(guess_xmid[peak_range])[iqr_index]         #interquartile range 
        
        binsize = 2.0*iqr*(2.0*n_sample)**(-1.0/3.0)        # appropriate bin size 
        print 'Freedman-Diaconis binsize = ', binsize
       
        # recompute histogram using appropriate binsize 
        n_bins = int((x_max-x_min)/binsize) 
        vlos_hist, v_binedges = np.histogram(combined_vlos, bins=n_bins, range=[x_min, x_max])
        xlow = v_binedges[:-1]
        xhigh = v_binedges[1:] 
        xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])

        peak_xrange = (xmid > peak_min) & (xmid < peak_max) 

        # amplitude of dLOS distribution 
        vlos_amp = np.mean(vlos_hist[(xmid >= -100.0) & (xmid < 100.0)])
        print vlos_amp, np.max(vlos_hist)
        
        p0 = [vlos_amp, 500.0]            # initial guess

        fa = {'x': xmid[peak_xrange], 'y': vlos_hist[peak_xrange]}
        if fit.lower() == 'expon': 
            peak_pars = mpfit.mpfit(mpfit_peak_expon, p0, functkw=fa, nprint=0)
        elif fit.lower() == 'gauss': 
            peak_pars = mpfit.mpfit(mpfit_peak_gauss, p0, functkw=fa, nprint=0)
        elif fit.lower() == 'loren': 
            peak_pars = mpfit.mpfit(mpfit_peak_loren, p0, functkw=fa, nprint=0)
        else: 
            raise NameError("Fit not yet coded") 
        
        bestfit_amp = peak_pars.params[0]
        sigma = peak_pars.params[1]

        # compute fpeak 
        xrange = (xmid >= -5000.0) & (xmid < 5000.0) 
        if fit.lower() == 'expon': 
            fpeak = np.sum(peak_expon(xmid[xrange], peak_pars.params))/np.sum(vlos_hist) 
        elif fit.lower() == 'gauss': 
            fpeak = np.sum(peak_gauss(xmid[xrange], peak_pars.params))/np.sum(vlos_hist) 
        elif fit.lower() == 'loren': 
            fpeak = np.sum(peak_loren(xmid[xrange], peak_pars.params))/np.sum(vlos_hist) 
        else: 
            raise NameError("Fit not yet coded") 

        print peak_expon(xmid[xrange], peak_pars.params)
        print np.sum(vlos_hist)

        # compute the chi-squared value for peak range
        if fit.lower() == 'expon': 
            E_peak = peak_expon(xmid, peak_pars.params)
        elif fit.lower() == 'gauss': 
            E_peak = peak_gauss(xmid, peak_pars.params)
        elif fit.lower() == 'loren': 
            E_peak = peak_gauss(xmid, peak_pars.params)
        else: 
            raise NameError("Fit not yet coded") 

        #print vlos_hist[peak_xrange]
        #print E_peak[peak_xrange]
        #print (vlos_hist[peak_xrange]-E_peak[peak_xrange])**2/E_peak[peak_xrange]
        chi2 = 1.0/np.float(len(xmid[peak_xrange])-3)*(np.sum((vlos_hist[peak_xrange]-E_peak[peak_xrange])**2 / E_peak[peak_xrange]))
        print 'Reduced Chi-Squared = ', chi2 

        if sanitycheck == True: 
            prettyplot() 
            pretty_colors = prettycolors() 
            sub = fig.add_subplot(3, 1, i_peak) 

            # fitting labels 
            if fit.lower() == 'expon': 
                fit_label = "Exponential "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
            elif fit.lower() == 'gauss': 
                fit_label = "Gaussian "+r"$(\sigma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
            elif fit.lower() == 'loren': 
                fit_label = "Lorentzian "+r"$(\gamma = "+str(peak_pars.params[1])+", A="+str(peak_pars.params[0])+")$"
        
            guess_scale = binsize/10.0
            sub.plot(guess_xmid, guess_scale*guess_vlos_hist, lw=2, color=pretty_colors[-1], label=r"$v_{LOS}$ binsize $=10.0$") 
            sub.plot(xmid, vlos_hist, lw=4, color=pretty_colors[0], label=r"$v_{LOS}$ binsize $="+str(binsize)+"$") 
            if fit.lower() == 'expon': 
                sub.plot(xmid, peak_expon(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
            elif fit.lower() == 'gauss': 
                sub.plot(xmid, peak_gauss(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
            elif fit.lower() == 'loren': 
                sub.plot(xmid, peak_loren(xmid, peak_pars.params), lw=4, color=pretty_colors[2], label=fit_label)
            else: 
                raise NameError("Fit not yet coded") 
        
            sub.text(-4750.0, 0.5*np.max(vlos_hist), 'peak: '+str(peak_min)+' - '+str(peak_max))
            sub.text(-4750.0, 0.43*np.max(vlos_hist), 'reduced chi2: '+str(chi2)) 
            if fit.lower() != 'loren': 
                sub.text(-1.0*sigma, 0.25*np.max(vlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 
            else: 
                sub.text(-1000.0, 0.25*np.max(vlos_hist), r"$f_{peak} = "+str(fpeak)+"$") 

            sub.set_xlabel(r"$v_{LOS}$", fontsize=20) 
            sub.set_xlim([-5000.0, 5000.0])
            sub.set_ylim([0.0, 1.25*np.max(vlos_hist)])
            sub.legend(loc='upper right', scatterpoints=1, prop={'size':14}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+catalog['name'].lower()+'_'+str(n_mocks)+'mocks_combined_vlos_peakfit_'+fit.lower()+'_peakrange_test.png', bbox_inches="tight")
    fig.clear() 

    return [sigma, fpeak]

def plot_vlos_peak(**cat_corr): 
    '''
    plot peak of vLOS distribution (testing appropriate binsize) 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    prettyplot() 
    pretty_colors = prettycolors()

    LOSdisp = vlos(**cat_corr) 
    LOSdisp.get_vlos_hist(binsize=100.0, normed=True) 
    peak_xrange =  (LOSdisp.xmid >= -3000.0) & (LOSdisp.xmid < 3000.0)         # approximate peak xrange
    #popt, pcov = curve_fit(LOSdisp.expon, (LOSdisp.xmid)[peak_xrange], (LOSdisp.vlos_hist)[peak_xrange])
    #print popt[0]
    dlos_amp = np.mean(LOSdisp.vlos_hist[peak_xrange])
    p0 = [dlos_amp, 5.0]            # initial guess
    fa = {'x': LOSdisp.xmid[peak_xrange], 'y': LOSdisp.vlos_hist[peak_xrange]}
    peak_pars = mpfit.mpfit(mpfit_peak_expon, p0, functkw=fa, nprint=0)
    
    bestfit_amp = peak_pars.params[0]
    sigma = peak_pars.params[1]
    print sigma

    fig = plt.figure(1) 
    sub = fig.add_subplot(111) 

    sub.plot(LOSdisp.xmid, LOSdisp.vlos_hist, color=pretty_colors[1]) 
    print LOSdisp.vlos_hist
    sub.plot(LOSdisp.xmid, peak_expon(LOSdisp.xmid, peak_pars.params), color=pretty_colors[5])
    sub.set_xlim([-5000.0, 5000.0])  
    #sub.set_ylim([0.0, 0.1])
    plt.show() 

def plot_fcpaper_vlos(cat_corrs): 
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
            mocks = range(1,46) 
        elif catalog['name'].lower() == 'cmass': 
            mocks = ['']
        else: 
            raise NameError('error')
    
        n_mocks = 0
        for i_mock in mocks: 
            i_catalog = catalog.copy()
            if catalog['name'].lower() != 'lasdamasgeo': 
                i_catalog['n_mock'] = i_mock 
                i_cat_corr = {'catalog':i_catalog, 'correction': correction} 
            else: 
                i_catalog['n_mock'] = int(''.join(list(i_mock)[0:-1]))
                i_catalog['letter'] = list(i_mock)[-1]
                i_cat_corr = {'catalog':i_catalog, 'correction': {'name': 'true'}} 
            
            los_disp_i = vlos(**i_cat_corr)         # import DLOS values from each mock 
            # Combine dLOS values 
            try: 
                combined_vlos
            except NameError: 
                combined_vlos = los_disp_i.vlos
            else: 
                combined_vlos = np.concatenate([combined_vlos, los_disp_i.vlos]) 

            n_mocks = n_mocks+1

        # Create histogram for combined dLOS values  (binsize is just guessed)
        x_min = -80000.0
        x_max = 80000.0
        binsize = 30.0 
        n_bins = int((x_max-x_min)/binsize) 
    
        vlos_hist, mpc_binedges = np.histogram(combined_vlos, bins=n_bins, range=[x_min, x_max], normed=True)
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

        sub.plot(xmid, vlos_hist, lw=lwid, color=cat_color, label=cat_label) 

        del combined_vlos
    
    #sub.set_title(cat_label) 
    sub.set_xlabel(r"v$_{\rm{LOS}}$ (km/sec)", fontsize=24) 
    sub.set_xlim([-5000.0, 5000.0])
    sub.set_ylim([0.0, 0.0005])
    #sub.set_xscale('log') 
    #sub.set_yscale('log') 
    sub.legend(loc='upper left', scatterpoints=1, prop={'size':20}) 

    fig_dir = fc_util.get_fig_dir() 
    fig.savefig(fig_dir+'fcpaper_vlos_dist.png', bbox_inches="tight")
    fig.clear() 

def random_sampling_test(): 
    fit_func = lambda x: np.exp(-0.5*x**2)

    for i in range(0, 1000):
        rand1 = np.random.random(1) 
        rand2 = np.random.random(1) 

        rand2 = (-3.0+rand2*6.0)
        peakpofr = fit_func(rand2) 

        while peakpofr <= rand1: 
            rand1 = np.random.random(1) 
            rand2 = np.random.random(1) 

            rand2 = (-3.0+rand2*6.0)
            peakpofr = fit_func(rand2) 
    
        try:
            sampled 
        except NameError: 
            sampled = rand2
        else: 
            sampled = np.concatenate([sampled, rand2])
    
    x_min = -3.0
    x_max = 3.0
    binsize = 0.1
    n_bins = int((x_max-x_min)/binsize) 
    x_hist, x_binedges = np.histogram(sampled, bins=n_bins, range=[x_min, x_max], density=True)
    xlow  = x_binedges[:-1]
    xhigh = x_binedges[1:] 
    xmid = np.array([0.5*(xlow[i]+xhigh[i]) for i in range(len(xlow))])
    print np.sum(x_hist) 
    fig = plt.figure(1) 
    sub = fig.add_subplot(111)
    sub.plot(xmid, x_hist) 
    sub.plot(xmid, fit_func(xmid)/np.sqrt(2*np.pi))
    plt.show() 

if __name__=="__main__": 
    #qpm_combined_vlos_zreal_comp()
    #random_sampling_test()
    cat_corr = {'catalog': {'name': 'qpm', 'n_mock':1}, 'correction': {'name': 'upweight'}} 
    combined_vlos_fit(fit='gauss', sanitycheck=True, **cat_corr)
    combined_vlos_fit(fit='expon', sanitycheck=True, **cat_corr)
    #combined_vlos_fit(fit='gauss', sanitycheck=True, vlos_method='vpec', **cat_corr)
    #combined_vlos_fit(fit='expon', sanitycheck=True, vlos_method='vpec', **cat_corr)

    #combined_vlos_fit(fit='gauss', sanitycheck=True, vlos_method='cDz', **cat_corr)
    #combined_vlos_fit(fit='expon', sanitycheck=True, vlos_method='cDz', **cat_corr)
    '''
    #combined_vlos_dist(**cat_corr)
    cat_corr = {'catalog': {'name': 'lasdamasgeo'}, 'correction': {'name': 'upweight'}} 
    #combined_vlos_dist(**cat_corr)
    combined_vlos_fit(fit='gauss', sanitycheck=True, **cat_corr)

    #cat_corr = {'catalog': {'name': 'qpm', 'n_mock':1}, 'correction': {'name': 'upweight'}} 
    #plot_vlos_peak(**cat_corr) 
    correction = {'name': 'upweight'}
    cat_corrs = [{'catalog': {'name':catalog}, 'correction': correction} for catalog in ['cmass', 'lasdamasgeo', 'qpm', 'tilingmock']]
    #plot_fcpaper_vlos(cat_corrs)
    '''
