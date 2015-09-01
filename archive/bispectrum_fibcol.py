import numpy as np
import pylab as py
import math as m
import matplotlib.pyplot as plt
from matplotlib import rc

'''
Analyze and plot power spectrum and bispectrum with fibercollision corrections!
'''
# Classes ---------------------------------------------
class bispec: 
    def __init__(self, mock, fibcolcorr): 
        ''' 
        Specify which mock catalog (tilingmock) 
        and fibercollision correction method  (true, delta, or peak) 
        '''
        if mock == 'tilingmock': 
            file_dir = '/mount/riachuelo1/hahn/bispec/tiling_mock/'
            file_prefix = 'BISP_cmass-boss5003sector-icoll012.'
            file_suffix = '.grid360.nmax.nstep3.P020000.box4000'
            
            # correction flags 
            if fibcolcorr == 'peak': 
                file_corr = 'fibcoll.dat.peak.sigma5.98.fpeak1.0'
            elif fibcolcorr == 'delta': 
                file_corr = 'fibcoll.dat.delta'
            elif fibcolcorr == 'true': 
                file_corr = 'dat'

            self.scale = 4000   # Mpc 
        self.file = ''.join([file_dir, file_prefix, file_corr, file_suffix]) 

    def readfile(self): 
        ''' 
        Read in bispectrum file given file name from init
        '''
        file_data = np.loadtxt(self.file) # load text 
        k_fund = (2.0*m.pi)/np.float(self.scale) # k_fundamental 
        self.kfund = k_fund
        self.k1 = k_fund*file_data[:,0]     # k1,k2,k3 of triangle 
        self.k2 = k_fund*file_data[:,1]
        self.k3 = k_fund*file_data[:,2]
        self.avgk = np.array([np.mean([self.k1[i], self.k2[i], self.k3[i]]) for i in range(len(self.k1))]) # mean(k1,k2,k3)
        self.Pk1 = file_data[:,3]           # power spectrum: P(k1), P(k2), P(k3) 
        self.Pk2 = file_data[:,4]
        self.Pk3 = file_data[:,5]
        self.Bk = file_data[:,6]            # Bispectrum B(k)
        self.Q = file_data[:,7]             # Reduced Bispectrum Q(k)

class powerspec: 
    def __init__(self, mock, fibcolcorr): 
        ''' 
        Specify which mock catalog (tilingmock) 
        and fibercollision correction method  (true, delta, or peak) 
        '''
        if mock == 'tilingmock': 
            file_dir = '/mount/riachuelo1/hahn/power/tiling_mocks/'
            file_prefix = 'power_cmass-boss5003sector-icoll012.'
            file_suffix = '.grid360.P020000.box4000'
            
            # correction flags
            if fibcolcorr == 'peak': 
                file_corr = 'fibcoll.dat.peak.sigma5.98.fpeak1.0'
            elif fibcolcorr == 'delta': 
                file_corr = 'fibcoll.dat.delta'
            elif fibcolcorr == 'true': 
                file_corr = 'dat'

            self.scale = 4000   # Mpc
        self.file = ''.join([file_dir, file_prefix, file_corr, file_suffix]) # filename 

    def readfile(self): 
        '''
        Read in P(k) file
        '''
        file_data = np.loadtxt(self.file) 
        self.k  = file_data[:,0]
        self.Pk = file_data[:,1]

# Functions ----------------------------------------
def residual(arr1, arr2): 
    '''
    Given two arrays of same dimensions, function returns array1/array2
    '''
    if len(arr1) != len(arr2): 
        raise TypeError("Input array lengths do not match.")
    else: 
        resid = np.array([arr1[i]/arr2[i] for i in range(len(arr1))])
    return resid

def residual_outlier_triangles(resid='Bk', correction='peak'): 
    '''
    obtain the triangle sides for B(k) or Q(k) residual outliers
    '''
    bisp_true = bispec('tilingmock', '')
    bisp_true.readfile()
    bisp_corr = bispec('tilingmock', correction)
    bisp_corr.readfile()
    if (resid=='Bk'): 
        corr_resid = residual(bisp_corr.Bk, bisp_true.Bk)
    elif (resid=='Qk'): 
        corr_resid = residual(bisp_corr.Q, bisp_true.Q)
    index_resid_outlier = (corr_resid > 1.5) | (corr_resid < 0.5)
    k_values = [bisp_corr.k1[index_resid_outlier]/bisp_corr.kfund, \
            bisp_corr.k2[index_resid_outlier]/bisp_corr.kfund, \
            bisp_corr.k3[index_resid_outlier]/bisp_corr.kfund]
    return k_values

def write_residual_outlier_triangles(k_values, resid='Bk', correction='peak'): 
    '''
    Writes out the triangle sides for kvalues output from residual_outlier_triangles
    '''
    file_dir = '/mount/riachuelo1/hahn/bispec/tiling_mock/'
    file_name = ''.join(['tilingmock_', resid, '_outlier_fibcol_', correction, '.dat'])
    f = open(''.join([file_dir, file_name]), 'w')
    print 'Writing ', ''.join([file_dir, file_name])
    for i in range(len(k_values[0])): 
        line = '\t'.join([str(k_values[0][i]), str(k_values[1][i]), str(k_values[2][i])])
        line = ''.join([line, '\n'])
        f.write(line)
    f.close()

# Plotting --------------------------------------------
def plot_Pk_resid(catalog='tilingmock'): 
    '''
    Plots P(k) residuals for mock catalog Fibercollision methods
    '''
    power_true = powerspec(catalog, 'true')
    #power_peak = powerspec(catalog, 'peak')
    power_delt = powerspec(catalog, 'delta')
    power_true.readfile()
    #power_peak.readfile()
    power_delt.readfile()
    
    bisp_true = bispec('tilingmock', 'true') 
    #bisp_peak = bispec('tilingmock', 'peak')
    bisp_delt = bispec('tilingmock', 'delta')
    bisp_true.readfile()
    #bisp_peak.readfile()
    bisp_delt.readfile()
    
    prettyplot()    # make plot pretty
    fig = plt.figure(1, figsize=(7,8))
    sub = fig.add_subplot(111)
    sub.plot(power_true.k, residual(power_delt.Pk, power_true.Pk), color='blue', 
            label=r'$\mathtt{P_{delta}(k)/P_{true}(k)}$')
    #sub.plot(power_true.k, residual(power_peak.Pk, power_true.Pk), color='red', 
    #        label=r'$\mathtt{P_{peak}(k)/P_{true}(k)}$')
    sub.scatter(bisp_true.k1, residual(bisp_delt.Pk1, bisp_true.Pk1), s=6, color='blue', 
            label=r'$\mathtt{P_{delta}(k1)/P_{true}(k1)}$')
    #sub.scatter(bisp_true.k1, residual(bisp_peak.Pk1, bisp_true.Pk1), s=6, color='red', 
    #        label=r'$\mathtt{P_{peak}(k1)/P_{true}(k1)}$')
    # specify axis 
    sub.set_xlim([10**-3,10**0])
    sub.set_ylim([0.5,1.5])
    sub.set_xlabel(r'k', fontsize=18)
    sub.set_ylabel(r'$\mathtt{P(k)}$', fontsize=18)
    sub.set_xscale('log')
    sub.grid(True)

    figdir = '/home/users/hahn/research/figures/boss/fiber_collision/tiling_mocks/'
    subleg = sub.legend(fontsize=18, loc='lower left', scatterpoints=1) # fig legend
    fig_name = ''.join([figdir, 'powerspec_residual_', catalog, '.png'])
    fig.savefig(fig_name)
    fig.clear()

def plot_Bk_fibcoll(catalog='tilingmock', triangle='all', kaxis='kavg', yaxis='Bk', resid='False'):
    '''
    Plot B(k) for Mock Catalog each Fibercollision methods with lots of options: 
    triangle: all/equilateral/
    kaxis: kavg or kmax
    yaxis: Bk or Bkavg
    resid: True or False   ; where you want to compute residuals or not 
    '''
    if resid == 'True':         # residuals 
        corr_label = [r'$B_{true}(k)$', r"$B_{delta}(k)/B_{true}(k)$"]#, r'$B_{peak}(k)/B_{true}(k)$']
        y_label = r"B"r"$(k)/$"r"B"r"$(k)_{true}$"
        ylimit = [0.5, 2.2]     # different yrange if residual
    else:                       # Not residuals 
        corr_label = [r'$B_{true}(k)$', r'$B_{delta}(k)$']#, r'$B_{peak}(k)$']
        y_label = r"B"r"$(k)$"
    
    fig = plt.figure(1, figsize=(7,8))          #figure
    sub = fig.add_subplot(111)

    corrections = ['true', 'delta']#, 'peak']     # fiber collision correction methods 
    corr_color = ['black', 'blue']#, 'red']
    for i_corr, corr in enumerate(corrections): 
        bisp = bispec(catalog, corr)            # read bispectrum
        bisp.readfile()

        if triangle == 'all':                   # for all triangles 
            # Compute k axis 
            if kaxis == 'kavg':                 # average(k1,k2,k3)
                k_values = [np.mean([bisp.k1[i], bisp.k2[i], bisp.k3[i]]) for i in range(len(bisp.k1))]
                k_label = r"$\mathtt{k_{avg}}$"
            elif kaxis == 'kmax':               # max(k1,k2,k3)
                k_values = [np.max([bisp.k1[i], bisp.k2[i], bisp.k3[i]]) for i in range(len(bisp.k1))]
                k_label = r"$\mathtt{k_{max}}$"
            Bk = bisp.Bk                        # B(k) 
        elif triangle == 'equilateral':         # Only equilateral triangles
            equalsides = (bisp.k1 == bisp.k2) & (bisp.k2 == bisp.k3)
            k_values = bisp.k1[equalsides]
            Bk = bisp.Bk[equalsides]
            if kaxis == 'kavg': 
                k_label = r"$\mathtt{k_{equal, avg}}$"
            elif kaxis == 'kmax': 
                k_label = r"$\mathtt{k_{equal, max}}$"

        if yaxis == 'Bkavg':                    # averages B(k) with same k value  
            k_unique = np.unique(k_values)
            Bk_avg = [np.mean(Bk[k_values == kval]) for kval in k_unique]
            k_values = k_unique
            Bk = Bk_avg
        elif yaxis == 'Bk':                     # all B(k) values 
            pass    # nothing 

        if resid == 'True':                     # B(k) residuals
            if corr == 'true': 
                Bk_true = Bk                    # no weights
            else: 
                sub.scatter(k_values, residual(Bk, Bk_true), s=8, color=corr_color[i_corr], label=corr_label[i_corr])
                resid_flag = '_resid'
        else:                                   # just B(k) 
            sub.scatter(k_values, Bk, s=8, color=corr_color[i_corr], label=corr_label[i_corr])
            resid_flag = ''
    # set axes
    sub.set_xlabel(k_label, fontsize=20)
    sub.set_ylabel(y_label, fontsize=20)
    sub.set_xlim([10**-3,10**0])
    sub.set_xscale('log')
    if catalog == 'tilingmock':
        catalog_display = 'Tiling Mock'
    if resid=='True': 
        sub.set_ylim(ylimit)
        sub.text(2*10**-3, ylimit[1]-0.25, catalog_display, fontsize=20)       # display catalog 
    subleg = sub.legend(fontsize=15, loc='lower left')
    
    fig_dir = '/home/users/hahn/research/figures/boss/fiber_collision/tiling_mocks/'    # figure directory
    fig_name = ''.join([fig_dir, 
        'bispec_', catalog, '_', triangle, 'triangles_', yaxis, resid_flag, '_', kaxis,'.png'])  # specify figure  
    print fig_name
    fig.savefig(fig_name)
    fig.clear()

def plot_Q123_fibcoll(catalog='tilingmock', triangle='all', kaxis='kavg', resid='True'): 
    '''
    Plots the reduced bispectrum Q_123(k) for Mocks Catalogs using different fibercollision methods 
    '''
    prettyplot()                                # pretty plot 
    fig = plt.figure(1, figsize=(7,8))          # configure figure 
    sub = fig.add_subplot(111)
    
    corrections = ['true', 'delta']#, 'peak']     # fiber collision correction methods 
    corr_color = ['black', 'blue', 'red']
    for i_corr, corr in enumerate(corrections): 
        bisp = bispec(catalog, corr)            # read bispectrum
        bisp.readfile()
        
        if triangle == 'all':                   # for all triangles 
            # Compute k axis 
            if kaxis == 'kavg':                 # average(k1,k2,k3)
                k_values =  bisp.avgk
                k_label = r"$\mathtt{k_{avg}}$"
            elif kaxis == 'kmax':               # max(k1,k2,k3)
                k_values = [np.max([bisp.k1[i], bisp.k2[i], bisp.k3[i]]) for i in range(len(bisp.k1))]
                k_label = r"$\mathtt{k_{max}}$"
            Qk = bisp.Q                        # Q(k) 
        elif triangle == 'equilateral':         # Only equilateral triangles
            equalsides = (bisp.k1 == bisp.k2) & (bisp.k2 == bisp.k3)
            k_values = bisp.k1[equalsides]
            Qk = bisp.Q[equalsides]
            if kaxis == 'kavg': 
                k_label = r"$\mathtt{k_{equal, avg}}$"
            elif kaxis == 'kmax': 
                k_label = r"$\mathtt{k_{equal, max}}$"

        corr_label = ''.join([r"$\mathtt{Q_{123, ", corr, "}(k)}$"])

        if resid == 'True':                     # B(k) residuals
            if corr == 'true': 
                Qk_true = Qk                    # no weights
            else: 
                corr_label = ''.join([r"$\mathtt{Q_{123,", corr, "}(k)/Q_{123, true}(k)}$"])
                sub.scatter(k_values, residual(Qk, Qk_true), s=8, color=corr_color[i_corr], label=corr_label)
                y_label = r"$\mathtt{Q_{123}(k)/Q_{123}(k)_{true}}$"
                ylimit = [0.5, 2.2]     # different yrange if residual
                resid_flag = '_resid'
        else:                                   # just B(k) 
            sub.scatter(k_values, Qk, s=8, color=corr_color[i_corr], label=corr_label)
            resid_flag = ''
            y_label = r"$\mathtt{Q_{123}(k)}$"
    # set axes
    sub.set_xlabel(k_label, fontsize=20)
    sub.set_ylabel(y_label, fontsize=20)
    sub.set_xlim([10**-3,10**0])
    sub.set_xscale('log')
    if catalog == 'tilingmock':
        catalog_display = 'Tiling Mock'                                     # set catalog text
    if resid=='True': 
        sub.set_ylim(ylimit)
        sub.text(2*10**-3, ylimit[1]-0.25, catalog_display, fontsize=20)       # display catalog text
    sub.legend(fontsize=15, loc='lower left')
    
    fig_dir = '/home/users/hahn/research/figures/boss/fiber_collision/tiling_mocks/'    # figure directory
    fig_filename = ''.join([fig_dir, 
        'bispec_', catalog, '_', triangle, 'triangles_Q123', resid_flag, '_', kaxis,'.png']) # specify figure 
    print fig_filename
    fig.savefig(fig_filename)
    fig.clear()

if __name__=="__main__":
    #plot_Pk_resid(catalog='tilingmock') 
    
    plot_Bk_fibcoll(triangle='all', kaxis='kavg', yaxis='Bk', resid='False')
    plot_Bk_fibcoll(triangle='all', kaxis='kavg', yaxis='Bkavg', resid='False')
    plot_Bk_fibcoll(triangle='all', kaxis='kavg', yaxis='Bk', resid='True')
    plot_Bk_fibcoll(triangle='all', kaxis='kavg', yaxis='Bkavg', resid='True')
    plot_Bk_fibcoll(triangle='equilateral', kaxis='kavg', yaxis='Bk', resid='True')
    plot_Bk_fibcoll(triangle='equilateral', kaxis='kavg', yaxis='Bkavg', resid='True')
    plot_Bk_fibcoll(triangle='all', kaxis='kmax', yaxis='Bk', resid='True')
    plot_Bk_fibcoll(triangle='all', kaxis='kmax', yaxis='Bkavg', resid='True')
    plot_Bk_fibcoll(triangle='equilateral', kaxis='kmax', yaxis='Bk', resid='True')
    plot_Bk_fibcoll(triangle='equilateral', kaxis='kmax', yaxis='Bkavg', resid='True')
    
    '''
    plot_Q123_fibcoll(catalog='tilingmock', triangle='all', kaxis='kavg', resid='False')
    plot_Q123_fibcoll(catalog='tilingmock', triangle='all', kaxis='kavg', resid='True')
    plot_Q123_fibcoll(catalog='tilingmock', triangle='all', kaxis='kmax', resid='False')
    plot_Q123_fibcoll(catalog='tilingmock', triangle='all', kaxis='kmax', resid='True')
    plot_Q123_fibcoll(catalog='tilingmock', triangle='equilateral', kaxis='kavg', resid='False')
    plot_Q123_fibcoll(catalog='tilingmock', triangle='equilateral', kaxis='kavg', resid='True')
    '''
