import numpy as np
import pylab as py
import math as m
import matplotlib.pyplot as plt
from matplotlib import rc

class bispec: 
    def __init__(self, mock, fibcolcorr): 
        if mock == 'tilingmock': 
            file_dir = '/mount/riachuelo1/hahn/bispec/tiling_mock/' 
            file_prefix = 'BISP_cmass-boss5003sector-icoll012.'
            file_suffix = '.grid360.nmax.nstep3.P020000.box4000'
            if fibcolcorr == 'peak': 
                file_corr = 'peaknbar.sigma5.978.fpeak1.0.fibcoll.dat.peaknbarcorr'
            elif fibcolcorr == 'delta': 
                file_corr = 'fibcoll.dat.deltacorr'
            else: 
                file_corr = 'dat'
            self.scale = 4000
        self.file = ''.join([file_dir, file_prefix, file_corr, file_suffix]) 
    def readfile(self): 
        file = np.loadtxt(self.file) 
        k_fund = (2.0*m.pi)/np.float(self.scale)
        self.kfund = k_fund
        self.k1 = k_fund*file[:,0]
        self.k2 = k_fund*file[:,1]
        self.k3 = k_fund*file[:,2]
        self.Pk1 = file[:,3]
        self.Pk2 = file[:,4]
        self.Pk3 = file[:,5]
        self.Bk = file[:,6]
        self.Q = file[:,7]
        self.avgk = np.array([np.mean([self.k1[i], self.k2[i], self.k3[i]]) for i in range(len(self.k1))])

class powerspec: 
    def __init__(self, mock, fibcolcorr): 
        if mock == 'tilingmock': 
            file_dir = '/mount/riachuelo1/hahn/power/tiling_mocks/'
            file_prefix = 'power_cmass-boss5003sector-icoll012.'
            file_suffix = '.grid360.P020000.box4000'
            if fibcolcorr == 'peak': 
                file_corr = 'fibcoll.dat.peaknbarcorr'
            elif fibcolcorr == 'delta': 
                file_corr = 'fibcoll.dat.deltacorr'
            else: 
                file_corr = 'dat'
            self.scale = 4000

        self.file = ''.join([file_dir, file_prefix, file_corr, file_suffix]) 
    def readfile(self): 
        file    = np.loadtxt(self.file) 
        self.k  = file[:,0]
        self.Pk = file[:,1]

def residual(arr1, arr2): 
    if len(arr1) != len(arr2): 
        raise TypeError("Input array lengths do not match.")
    else: 
        resid = np.array([arr1[i]/arr2[i] for i in range(len(arr1))])
    return resid

def plot_tilingmock_Pk_resid(): 
    '''
    Plots P(k) residuals for Tiling Mock Fibercollision methods
    '''
    figdir = '/home/users/hahn/research/figures/boss/bispectrum/fibcol/'
    power_true = powerspec('tilingmock', '')
    power_peak = powerspec('tilingmock', 'peak')
    power_delt = powerspec('tilingmock', 'delta')
    power_true.readfile()
    power_peak.readfile()
    power_delt.readfile()
    
    bisp_true = bispec('tilingmock', '') 
    bisp_peak = bispec('tilingmock', 'peak')
    bisp_delt = bispec('tilingmock', 'delta')
    bisp_true.readfile()
    bisp_peak.readfile()
    bisp_delt.readfile()

    fig5 = plt.figure(5, figsize=(7,8))
    sub5 = fig5.add_subplot(111)
    sub5.plot(power_true.k, residual(power_delt.Pk, power_true.Pk), \
            color='black', label=r'$P_{delta}(k)/P_{true}(k)$')
    sub5.plot(power_true.k, residual(power_peak.Pk, power_true.Pk), \
            color='red', label=r'$P_{peak}(k)/P_{true}(k)$')
    sub5.scatter(bisp_true.k1, residual(bisp_delt.Pk1, bisp_true.Pk1), \
            s=6, color='black', label=r'$P_{delta}(k1)/P_{true}(k1)$')
    sub5.scatter(bisp_true.k1, residual(bisp_peak.Pk1, bisp_true.Pk1), \
            s=6, color='red', label=r'$P_{peak}(k1)/P_{true}(k1)$')
    sub5.set_xlim([10**-3,10**0])
    sub5.set_ylim([0.5,1.5])
    sub5.set_xlabel(r'$k$', fontsize=15)
    sub5.set_ylabel(r'$P(k)$', fontsize=15)
    sub5.set_xscale('log')
    sub5.grid(True)
    sub5leg = sub5.legend(fontsize=15, loc='lower left', fancybox=True)
    sub5leg.get_frame().set_alpha(0.5)
    fig_name = ''.join([figdir, 'powerspec_tilingmock_pk_residual.png'])
    fig5.savefig(fig_name)

def plot_tilingmock_Bk(triangle='none', kaxis='kavg', yaxis='Bk', resid='False'):
    '''
    Plots B(k) for Tiling Mock Fibercollision methods, also computes residuals
    '''
    figdir = '/home/users/hahn/research/figures/boss/bispectrum/fibcol/'
    corrections = ['', 'delta', 'peak'] 
    corr_color = ['black', 'blue', 'red']
    if resid=='False': 
        corr_label = [r'$B_{true}(k)$', r'$B_{delta}(k)$', r'$B_{peak}(k)$']
        y_label = r"B"r"$(k)$"
    else: 
        corr_label = [r'$B_{true}(k)$', r"$B_{delta}(k)/B_{true}(k)$", r'$B_{peak}(k)/B_{true}(k)$']
        y_label = r"B"r"$(k)/$"r"B"r"$(k)_{true}$"
        ylimit = [0.5, 2.2]
    for i_corr, corr in enumerate(corrections): 
        # Load bispectrum
        bisp = bispec('tilingmock', corr) 
        bisp.readfile()
        # Load figure 
        fig4 = plt.figure(4, figsize=(7,8))
        if triangle=='none': 
            sub4 = fig4.add_subplot(111)
            # Compute k axis 
            if (kaxis=='kavg'): # average(k1,k2,k3)
                k_axis = [np.mean([bisp.k1[i], bisp.k2[i], bisp.k3[i]]) \
                        for i in range(len(bisp.k1))]
                k_label = r"k"r"$_{avg}$"
            elif (kaxis=='kmax'): # max(k1,k2,k3)
                k_axis = [np.max([bisp.k1[i], bisp.k2[i], bisp.k3[i]]) \
                        for i in range(len(bisp.k1))]
                k_label = r"k"r"$_{max}$"
            Bk = bisp.Bk
            kflag = kaxis
        elif triangle=='equilateral': # Only equilateral triangles
            equalsides = (bisp.k1 == bisp.k2) & (bisp.k2 == bisp.k3)
            k_axis = bisp.k1[equalsides]
            Bk = bisp.Bk[equalsides]
            k_label = r"k"r"$_{equal}$"
            kflag = 'equilateral'
        if (yaxis=='Bkavg'):
            k_unique = np.unique(k_axis)
            Bk_avg = []
            for kval in k_unique: 
                Bk_avg.append(np.mean(Bk[k_axis==kval]))
            k_axis = k_unique
            Bk = Bk_avg
        yaxisflag = yaxis 
        if (resid=='False'):
            sub4.scatter(k_axis, Bk, s=8, color=corr_color[i_corr], label=corr_label[i_corr])
            resid_flag = ''
        else: 
            if (corr==''): 
                Bk_true = Bk
            else: 
                sub4.scatter(k_axis, residual(Bk, Bk_true), s=8, color=corr_color[i_corr], label=corr_label[i_corr])
                resid_flag = '_resid'
    sub4.set_xlabel(k_label, fontsize=20)
    sub4.set_ylabel(y_label, fontsize=20)
    sub4.set_xlim([10**-3,10**0])
    if resid=='True': 
        sub4.set_ylim(ylimit)
        sub4.text(2*10**-3, ylimit[1]-0.25, "Tiling Mock", fontsize=20)
    sub4.set_xscale('log')
    sub4leg = sub4.legend(fontsize=15, loc='lower left', fancybox=True)
    sub4leg.get_frame().set_alpha(0.5)
    fig_name = ''.join([figdir, 'bispec_tilingmock_', yaxisflag, resid_flag, '_', kflag,'.png'])
    print fig_name
    fig4.savefig(fig_name)
    fig4.clear()

def plot_tilingmock_Q123_resid(triangle='none'): 
    '''
    Plots Q123 for Tiling Mocks using different fibercollision methods 
    '''
    figdir = '/home/users/hahn/research/figures/boss/bispectrum/fibcol/'
    fig5 = plt.figure(5, figsize=(7,8))
    sub = fig5.add_subplot(111)
    # Load correction method 
    correction = ['delta', 'peak']
    corrcolor = ['black', 'red']
    bisp_true = bispec('tilingmock', '')
    bisp_true.readfile()
    for i_corr, method in enumerate(correction):
        bisp_corr = bispec('tilingmock', method)
        bisp_corr.readfile()
        # for different specified triangles 
        if (triangle=='none'): 
            sub.scatter(bisp_true.avgk, residual(bisp_corr.Q, bisp_true.Q),\
                    s=4, color=corrcolor[i_corr], label=r''.join(['$B_{', method, '}(k)/B_{true}(k)$']))
            #print np.max(residual(bisp_corr.Q, bisp_true.Q)), np.min(residual(bisp_corr.Q, bisp_true.Q))
            figname = 'bispec_tilingmock_Qk_residual.png'
        if (triangle=='equilateral'): 
            equalsides = (bisp_true.k1 == bisp_true.k2) & (bisp_true.k2 == bisp_true.k3)
            equalk = bisp_true.k1[equalsides]
            sub.scatter(equalk, residual(bisp_corr.Q[equalsides], bisp_true.Q[equalsides]), \
                    s=4, color=corrcolor[i_corr], label=r''.join(['$B_{', method, '}(k)/B_{true}(k)$']))
            #print np.max(residual(bisp_corr.Q[equalsides], bisp_true.Q[equalsides])), \
                    #        np.min(residual(bisp_corr.Q[equalsides], bisp_true.Q[equalsides]))
            figname = 'bispec_tilingmock_Qk_residual_equilateral.png'
            sub.text(2*10**-3, 1.65, "Equilateral Triangles Only", fontsize=15)
    sub.set_xlim([10**-3, 10**0])
    sub.set_ylim([0.0, 2.0])
    sub.set_xlabel(r"k"r"$_{avg}$", fontsize=20)
    sub.set_ylabel(r"$Q/Q_{true}$", fontsize=20)
    sub.set_xscale('log')
    sub.text(2*10**-3, 1.75, "Tiling Mock", fontsize=20)
    subleg = sub.legend(fontsize=15, loc='lower left', fancybox=True)
    fig_name = ''.join([figdir, figname])
    fig5.savefig(fig_name)
    fig5.clear()

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

if __name__=="__main__":
    rc('text', usetex=True)
    rc('font', family='serif')

#    peak_outlier_triangles = residual_outlier_triangles(resid='Bk', correction='peak')
#    delta_outlier_triangles = residual_outlier_triangles(resid='Bk', correction='delta')
#    print 'Number of outliers in Delta Correction =', len(delta_outlier_triangles[0])
#    print 'Number of outliers in Peak Correction =', len(peak_outlier_triangles[0])
#    write_residual_outlier_triangles(peak_outlier_triangles, resid='Bk', correction='peak')
#    write_residual_outlier_triangles(delta_outlier_triangles, resid='Bk', correction='delta')

    plot_tilingmock_Bk(triangle='none', kaxis='kavg')
    plot_tilingmock_Bk(triangle='none', kaxis='kmax')
    plot_tilingmock_Bk(triangle='none', kaxis='kavg', yaxis='Bkavg')
    plot_tilingmock_Bk(triangle='none', kaxis='kmax', yaxis='Bkavg')
    plot_tilingmock_Bk(triangle='none', kaxis='kavg', yaxis='Bk', resid='True')
    plot_tilingmock_Bk(triangle='none', kaxis='kmax', yaxis='Bk', resid='True')
    plot_tilingmock_Bk(triangle='none', kaxis='kavg', yaxis='Bkavg', resid='True')
    plot_tilingmock_Bk(triangle='none', kaxis='kmax', yaxis='Bkavg', resid='True')
#    plot_tilingmock_Bk_resid(triangle='none', plotflag='Bk')
#    plot_tilingmock_Bk_resid(triangle='none', plotflag='Bkamplitude')
#    plot_tilingmock_Bk_resid(triangle='equilateral', plotflag='Bk')
#    plot_tilingmock_Bk_resid(triangle='equilateral', plotflag='Bkamplitude')
#    plot_tilingmock_Q123_resid(triangle='none')
#    plot_tilingmock_Q123_resid(triangle='equilateral')
