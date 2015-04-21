import numpy as np
import pylab as py
from scipy.optimize import curve_fit
from matplotlib import rc
import sys

class dlos_histogram:
    def __init__(self): 
        pass
    def readFile(self, catalog='LasDamasGeo', version='1a', n_mock=1):
        '''
        Read in dLOS values for specified catalog and file
        Catalogs: LasDamasGeo, LasDamasBox, PTHalo
        '''
        if catalog=='LasDamasGeo': 
            dlos_file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
            lasdamas_n_mock = ''
            lasdamas_letter = ''
            for i in n_mock: 
                if i.isdigit(): 
                    lasdamas_n_mock += i 
                else: 
                    lasdamas_letter += i 
            dlos_file_name = 'DLOS_sdssmock_gamma_lrgFull_zm_oriana'+str(int(lasdamas_n_mock)+100)[1:3]+lasdamas_letter+'_no.rdcz.dat'
            version = ''
        elif catalog=='LasDamasBox': 
            #if len(input) !=1: 
            #    raise NameError('LasDamasBox needs [number]') 
            #if type(input[0]) is not int: 
            #    raise NameError('LasDamasBox needs [number]') 
            dlos_file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Box/'
        elif catalog == 'PTHalo': 
            if type(version) is not str: 
                raise NameError("PTHalo requires version_string (e.g. 'v11p0')") 
            if type(n_mock) is not int: 
                raise NameError("PTHalo requires mock file number (e.g. 1)") 
            
            dlos_file_dir = ''.join(['/mount/riachuelo1/hahn/data/PTHalo/', version, '/'])
            dlos_file_name = ''.join([
                'DLOS_cmass_dr11_north_ir4', str(n_mock+1000)[1:4], '.', 
                '.'.join(version.split('p')), '.wghtv.txt'
                ])
        elif catalog == 'TilingMock': 
            dlos_file_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'
            dlos_file_name = 'DLOS_cmass-boss5003sector-icoll012.dat'
            version = ''
            n_mock = ''
        elif catalog.lower() == 'qpm': 
            dlos_file_dir = '/mount/riachuelo1/hahn/data/QPM/dr12/ngc/'     # currently hardcoded for just DR12 NGC 
            dlos_file_name = 'DLOS_a0.6452_'+str(n_mock+10000)[1:5]+'.dr12c_ngc.rdz.dat'
            version = 'dr12'
        self.catalog = catalog  # catalog
        self.version = version   # File specifier
        self.n_mock = n_mock
        self.file = ''.join([dlos_file_dir, dlos_file_name])
        self.dlos = np.loadtxt(self.file)
    def get_hist(self): 
        ''' 
        Compute the histogram for specified dLOS file
        '''
        fig = py.figure(1)
        dump = fig.add_subplot(111)
        mpc_bin = -1000.0+0.1*np.array(range(20001))
        hist_dlos = dump.hist(self.dlos, mpc_bin)
        self.xlow = [hist_dlos[1][i] for i in range(len(hist_dlos[1])-1)] 
        self.xhigh = [hist_dlos[1][i+1] for i in range(len(hist_dlos[1])-1)] 
        self.xmid = [(hist_dlos[1][i]+hist_dlos[1][i+1])/2.0 for i in range(len(hist_dlos[1])-1)] 
        self.hist = hist_dlos[0]
        fig.clear()
    def expon(self, x, sig):
        return np.max(self.hist)*np.exp(-x/sig)
    def gauss(self, x, sig):
        return np.max(self.hist)*np.exp(-0.5*x**2/sig**2)

def plot_curvefit_dlos(catalog='PTHalo', version='v11p0', n_mock=1, fit='expon'): 
    '''
    plot the bestfit exponential/gaussian to dLOS histogram 
    dlos xrange is currently hardcoded
    '''
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

def get_dlos_curvefit_sigma(catalog='PTHalo', version='v11p0', n_mock=1, fit='expon'): 
    '''
    calculate sigma of the bestfit exponential/gaussian to dLOS histogram 
    dlos xrange is currently hardcoded
    '''
    dlos_hist = dlos_histogram()
    dlos_hist.readFile(catalog=catalog, version=version, n_mock=n_mock) 
    dlos_hist.get_hist()

    if fit=='expon': 
        popt, pcov = curve_fit(dlos_hist.expon, np.array(dlos_hist.xmid[10000:11000]), np.array(dlos_hist.hist[10000:11000]))
    elif fit=='gauss':
        popt, pcov = curve_fit(dlos_hist.gauss, np.array(dlos_hist.xmid[9500:10500]), np.array(dlos_hist.hist[9500:10500]))
    return popt[0]

def average_sigma(catalog='PTHalo', version='v11p0', n_mock=range(1,11), fit='expon'): 
    '''
    calculate average(sigma) of the bestfit exponential/gaussian of fiber-collided pair dLOS for mock catalogs
    '''
    sigma_sum = 0.0 
    for n_file in n_mock: 
        sigma_n = get_dlos_curvefit_sigma(catalog=catalog, version=version, n_mock=n_file, fit=fit) 
        print catalog, version, n_file, sigma_n
        sigma_sum = sigma_sum+sigma_n
    print sigma_sum 
    return sigma_n/float(len(n_mock)) 

if __name__=="__main__": 
    #plot_curvefit_dlos(catalog='LasDamasGeo', n_mock='1a', fit='expon')
    #plot_curvefit_dlos(catalog='PTHalo', version='v11p0', n_mock=1, fit='expon')
    #plot_curvefit_dlos(catalog='PTHalo', version='v11p0', n_mock=2, fit='expon')
    #plot_curvefit_dlos(catalog='TilingMock', fit='expon')
    plot_curvefit_dlos(catalog='QPM', fit='expon')
    #print "Average Sigma_exp for PTHALO 1-10", average_sigma(catalog='PTHalo', version='v11p0', n_mock=range(1,11), fit='expon')
