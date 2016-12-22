'''


Code for dealing fiber collisions in the bispectrum.


'''
import numpy as np 
from corr_spec.corr_spec import CorrSpec 
from corr_spec.corr_multi import build_multipro

# Plotting
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt
#from ChangTools.plotting import png2pdf 
from ChangTools.plotting import prettyplot 
from ChangTools.plotting import prettycolors 
from ChangTools.codif import notif

def MeasureBisp(catalog='qpm', correction='upweight', n_mocks=1, Ngrid=360, Nthreads=1): 
    ''' Wrapper for measuring the bispectrum
    '''
    if isinstance(n_mocks, list): 
        n_list = n_mocks
    else: 
        n_list = range(1, n_mocks+1)
    build_multipro('bk', catalog, correction, n_list, Nthreads=Nthreads, Ngrid=Ngrid)
    return None


def Bisp_FcImpact(catalog='nseries', GM_klim=False, quadrupole=False): 
    ''' Quantify the impact of fiber collisions on the bispectrum.
    Compare the bispectrum of the upweighted catalog WITH fiber collisions 
    to the bispectrum of the true catalog.
    '''
    if catalog == 'nseries': 
        n_mocks = 84
    elif catalog == 'qpm': 
        n_mocks = 100
    else: 
        raise NotImplementedError
    
    for i_mock in range(1, n_mocks+1): 
        true_dict = { 
                'catalog': {'name': catalog, 'n_mock': i_mock}, 
                'correction': {'name': 'true'}
                }
        bisp_true = CorrSpec('bk', true_dict, Ngrid=360)
        bisp_true.read()
        
        upw_dict = { 
                'catalog': {'name': catalog, 'n_mock': i_mock}, 
                'correction': {'name': 'upweight'}
                }
        bisp_upw = CorrSpec('bk', upw_dict, Ngrid=360)
        bisp_upw.read()

        if i_mock == 1: 
            if not quadrupole: 
                sum_B_true = bisp_true.bk
                sum_B_upw = bisp_upw.bk
                sum_B_ratio = bisp_upw.bk/bisp_true.bk
            else: 
                sum_B_true = bisp_true.b2k
                sum_B_upw = bisp_upw.b2k
                sum_B_ratio = bisp_upw.b2k/bisp_true.b2k
            k1 = bisp_upw.k1
            k2 = bisp_upw.k2
            k3 = bisp_upw.k3
        else: 
            if not quadrupole: 
                sum_B_true += bisp_true.bk
                sum_B_upw += bisp_upw.bk
                sum_B_ratio += bisp_upw.bk/bisp_true.bk
            else: 
                sum_B_true += bisp_true.b2k
                sum_B_upw += bisp_upw.b2k
                sum_B_ratio += bisp_upw.b2k/bisp_true.b2k

    avg_B_true = sum_B_true / np.float(n_mocks)
    avg_B_upw = sum_B_upw / np.float(n_mocks) 
    avg_B_ratio = sum_B_ratio / np.float(n_mocks) 

    if GM_klim: 
        krange = np.where(
                (k1 >= 0.03) & (k1 <= 0.22) & 
                (k2 >= 0.03) & (k2 <= 0.22) & 
                (k3 >= 0.03) & (k3 <= 0.22) 
                )
        avg_B_true = avg_B_true[krange]
        avg_B_upw = avg_B_upw[krange]
        k1 = k1[krange]
        k2 = k2[krange]
        k3 = k3[krange]
        print k3.min(), k3.max() 
    
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(15,10))
    sub = fig.add_subplot(211)

    sub.scatter(range(len(avg_B_true)), avg_B_true, c=pretty_colors[0], lw=0, s=9, label='True')
    sub.scatter(range(len(avg_B_upw)), avg_B_upw, c=pretty_colors[3], lw=0, s=3, label='Upweighted')
    sub.set_xlim([-100, len(avg_B_true)+10])
    
    sub.set_xticklabels([])
    if not quadrupole: 
        sub.set_ylim([10.**7, 1.*10**10])
        sub.set_ylabel(r'$\mathtt{B(k)}$', fontsize=25) 
    else: 
        sub.set_ylim([10.**7, 1.*10**10])
        sub.set_ylabel(r'$\mathtt{B_2(k)}$', fontsize=25) 
    sub.set_yscale("log") 
    sub.legend(loc='upper right', scatterpoints=1) 

    sub = fig.add_subplot(212)
    
    sub.scatter(range(len(avg_B_true)), avg_B_upw/avg_B_true, c=pretty_colors[0], lw=0, s=7, label='Other')
    print 'All triangles'
    print min(avg_B_ratio), max(avg_B_ratio)
    print 'Minimum ratio ', min(avg_B_upw/avg_B_true)
    print 'Maximum ratio', max(avg_B_upw/avg_B_true) 

    for i_t, triangle in enumerate(['equilateral', 'acute', 'obtuse', 'extended']): 
        tri_index = ClassifyTriangles(k1, k2, k3, triangle=triangle)
        sub.scatter(tri_index[0], avg_B_upw[tri_index]/avg_B_true[tri_index], 
                c=pretty_colors[i_t+2], lw=0, s=7, label=str(triangle.title()))
        print triangle.title() 
        print min(avg_B_ratio[tri_index]), max(avg_B_ratio[tri_index])
        print 'Minimum ratio ', min(avg_B_upw[tri_index]/avg_B_true[tri_index])
        print 'Maximum ratio', max(avg_B_upw[tri_index]/avg_B_true[tri_index]) 

    sub.legend(loc='upper left', scatterpoints=1) 

    sub.set_xlim([-100, len(avg_B_true)+10])
    sub.set_xlabel('Triangle Index', fontsize=25) 
    if not quadrupole: 
        sub.set_ylim([0.975, 1.15])
    else: 
        sub.set_ylim([0.8, 1.5])
    sub.set_ylabel(r'$\mathtt{B^{upw}(k)/B^{true}(k)}$', fontsize=25) 

    if GM_klim: 
        gm_str = '.GM_krange'
    else: 
        gm_str = ''

    if quadrupole: 
        quad_str = '.quadrupole'
    else: 
        quad_str = ''

    fig_name = ''.join(['figure/',
        'Bispectrum_FcImpact', 
        '.', catalog, 
        gm_str, 
        quad_str, 
        '.png']) 
    fig.savefig(fig_name, bbox_inches='tight') 
    plt.close()
    #png2pdf(fig_name)
    return None


def Bisp_CMASScomp(catalog='nseries', GM_klim=False, quadrupole=False): 
    ''' Compare mock bispectrum to CMASS bispectrum
    '''
    if catalog == 'nseries': 
        n_mocks = 84
    elif catalog == 'qpm': 
        n_mocks = 100
    else: 
        raise NotImplementedError
    
    for i_mock in range(1, n_mocks+1): 
        upw_dict = { 
                'catalog': {'name': catalog, 'n_mock': i_mock}, 
                'correction': {'name': 'upweight'}
                }
        bisp_upw = CorrSpec('bk', upw_dict, Ngrid=360)
        bisp_upw.read()

        if i_mock == 1: 
            if not quadrupole: 
                sum_B_upw = bisp_upw.bk
            else: 
                sum_B_upw = bisp_upw.b2k
            k1 = bisp_upw.k1
            k2 = bisp_upw.k2
            k3 = bisp_upw.k3
        else: 
            if not quadrupole: 
                sum_B_upw += bisp_upw.bk
            else: 
                sum_B_upw += bisp_upw.b2k
    avg_B_upw = sum_B_upw / np.float(n_mocks) 

    # import CMASS bispectrum
    cmass_file = ''.join([
        '/mount/riachuelo1/hahn/power/CMASS/', 
        'BISPv5_cmass-dr12v4-N-Reid_fidcomso.dat.grid360.nmax40.ncut3.s3.P020000.box3600' 
        ]) 
    cmass_k_fund = (2.0*np.pi)/3600.        # k fundamental 
    k1_c, k2_c, k3_c, bk_c, b2k_c = np.loadtxt(cmass_file, unpack=True, usecols=[0, 1, 2, 6, 8])
    k1_c *= cmass_k_fund
    k2_c *= cmass_k_fund
    k3_c *= cmass_k_fund

    if not quadrupole: 
        B_cmass = bk_c
    else: 
        B_cmass = b2k_c

    if not np.array_equal(k1_c, k1): 
        raise ValueError

    if GM_klim: 
        krange = np.where(
                (k1 >= 0.03) & (k1 <= 0.22) & 
                (k2 >= 0.03) & (k2 <= 0.22) & 
                (k3 >= 0.03) & (k3 <= 0.22) 
                )
        avg_B_upw = avg_B_upw[krange]
        k1 = k1[krange]
        k2 = k2[krange]
        k3 = k3[krange]

        B_cmass = B_cmass[krange]
    
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(15,5))
    sub = fig.add_subplot(111)
    sub.scatter(range(len(avg_B_upw)), avg_B_upw, c=pretty_colors[3], lw=0, s=3, label=catalog.title())
    sub.scatter(range(len(B_cmass)), B_cmass, c=pretty_colors[0], lw=0, s=2, label='CMASS')
    sub.set_xlim([-100, len(avg_B_upw)+10])
    sub.set_xlabel('Triangle Index', fontsize=25)
    if not quadrupole: 
        sub.set_ylim([10.**7, 1.*10**10])
        sub.set_ylabel(r'$\mathtt{B(k)}$', fontsize=25) 
    else: 
        sub.set_ylim([10.**7, 1.*10**10])
        sub.set_ylabel(r'$\mathtt{B_2(k)}$', fontsize=25) 
    sub.set_yscale("log") 
    sub.legend(loc='upper right', scatterpoints=1) 
    
    if GM_klim: 
        gm_str = '.GM_krange'
    else: 
        gm_str = ''

    if quadrupole: 
        quad_str = '.quadrupole'
    else: 
        quad_str = ''

    fig_name = ''.join(['figure/',
        'Bispectrum_CMASScomp', 
        '.', catalog, 
        gm_str, 
        quad_str, 
        '.png']) 
    fig.savefig(fig_name, bbox_inches='tight') 
    plt.close()
    #png2pdf(fig_name)
    return None


def BispTriangle(catalog='nseries', GM_klim=False, quadrupole=False):
    ''' Quantify the impact of fiber collisions on the bispectrum.
    Compare the bispectrum of the upweighted catalog WITH fiber collisions 
    to the bispectrum of the true catalog.
    '''
    if catalog == 'nseries': 
        n_mocks = 10 
    elif catalog == 'qpm': 
        n_mocks = 100
    else: 
        raise NotImplementedError
    
    for i_mock in range(1, n_mocks+1): 
        true_dict = { 
                'catalog': {'name': catalog, 'n_mock': i_mock}, 
                'correction': {'name': 'true'}
                }
        bisp_true = CorrSpec('bk', true_dict, Ngrid=360)
        bisp_true.read()
        
        upw_dict = { 
                'catalog': {'name': catalog, 'n_mock': i_mock}, 
                'correction': {'name': 'upweight'}
                }
        bisp_upw = CorrSpec('bk', upw_dict, Ngrid=360)
        bisp_upw.read()

        if i_mock == 1: 
            if not quadrupole: 
                sum_B_true = bisp_true.bk
                sum_B_upw = bisp_upw.bk
                sum_B_ratio = bisp_upw.bk/bisp_true.bk
            else: 
                sum_B_true = bisp_true.b2k
                sum_B_upw = bisp_upw.b2k
                sum_B_ratio = bisp_upw.b2k/bisp_true.b2k
            k1 = bisp_upw.k1
            k2 = bisp_upw.k2
            k3 = bisp_upw.k3
        else: 
            if not quadrupole: 
                sum_B_true += bisp_true.bk
                sum_B_upw += bisp_upw.bk
                sum_B_ratio += bisp_upw.bk/bisp_true.bk
            else: 
                sum_B_true += bisp_true.b2k
                sum_B_upw += bisp_upw.b2k
                sum_B_ratio += bisp_upw.b2k/bisp_true.b2k

    avg_B_true = sum_B_true / np.float(n_mocks)
    avg_B_upw = sum_B_upw / np.float(n_mocks) 
    avg_B_ratio = sum_B_ratio / np.float(n_mocks) 

    k2k1 = k2/k1
    k3k1 = k3/k1
    x_bins = np.arange(k2k1.min(), k2k1.max()+0.02, 0.02)
    y_bins = np.arange(k3k1.min(), k3k1.max()+0.02, 0.02)
    true_bisp_grid = np.zeros((len(x_bins)-1, len(y_bins)-1))
    upw_bisp_grid = np.zeros((len(x_bins)-1, len(y_bins)-1))

    for i_x in range(len(x_bins)-1): 
        for i_y in range(len(y_bins)-1): 
            lim = np.where(
                    (k2k1 >= x_bins[i_x]) & (k2k1 < x_bins[i_x+1]) &
                    (k3k1 >= y_bins[i_y]) & (k3k1 < y_bins[i_y+1]))

            true_bisp_grid[i_x, i_y] = np.sum(avg_B_true[lim])
            upw_bisp_grid[i_x, i_y] = np.sum(avg_B_upw[lim])
    
    bisp_grid = upw_bisp_grid / true_bisp_grid
    #cm = plt.cm.get_cmap('RdYlBu')
    #sc = plt.scatter(k3k1, k2k1, c=avg_B_true, s=10, lw=0, vmin=0, vmax=1, cmap=cm)
    #plt.colorbar(sc)
    X, Y = np.meshgrid(x_bins[:-1], y_bins[:-1])
    #print np.max(np.log10(bisp_grid)[np.where(np.isnan(np.log10(bisp_grid)) == False)])
    #print np.min(np.log10(bisp_grid)[np.where((np.isnan(np.log10(bisp_grid)) == False) & (bisp_grid > 0.))])
    print np.max(bisp_grid[np.where(np.isnan(bisp_grid) == False)])
    print np.min(bisp_grid[np.where(np.isnan(bisp_grid) == False)])
    #print np.min(np.log10(bisp_grid)[np.where((np.isnan(np.log10(bisp_grid)) == False) & (bisp_grid > 0.))])
    #plt.pcolormesh(np.log10(bisp_grid), vmin=8.3, vmax=10.8)
    #plt.imshow(np.log10(bisp_grid), vmin=8.3, vmax=10.8, interpolation='bilinear')
    prettyplot() 
    pretty_colors = prettycolors() 
    fig = plt.figure(figsize=(10,3))
    sub = fig.add_subplot(121)
    sub.set_xlim([0.05, 1.05]) 
    sub.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0]) 
    sub.set_ylim([0.45, 1.05]) 
    sub.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.]) 
    sub.set_ylabel(r'$\mathtt{k_2/k_1}$', fontsize=24)
    sub = fig.add_subplot(122)
    sub.set_xlim([0.05, 1.05]) 
    sub.set_xticks([0.2, 0.4, 0.6, 0.8, 1.0]) 
    sub.set_ylim([0.45, 1.05]) 
    sub.set_yticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.]) 
    sub.set_yticklabels([])
    bkgd = fig.add_subplot(111, frameon=False)
    bkgd.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    bkgd.set_xlabel(r'$\mathtt{k_3/k_1}$', fontsize=24)
    
    fig.subplots_adjust(wspace=0.05, hspace=0.)
    #bisp_grid[np.where(bisp_grid == 0)] = -np.inf
    #bplot = plt.pcolormesh(bisp_grid, vmin=0.98, vmax=1.1, cmap=cm.OrRd)#, interpolation='bilinear')
    #cbar = plt.colorbar(bplot, orientation='horizontal') 
    #cbar.set_ticks([1., 1.05, 1.1]) # add the labels

    #plt.imshow(bisp_grid, vmin=., vmax=1.1, interpolation='bilinear')

    
    fig_name = ''.join(['figure/',
        'BispTriangle_FcImpact', '.png']) 
    plt.savefig(fig_name, bbox_inches='tight') 
    plt.close()
    #png2pdf(fig_name)
    return None


def ClassifyTriangles(k1, k2, k3, triangle='equilateral'):
    '''
    Given k1, k2, k3, returns indices for (k1,k2,k3) that satify the specified triangle type
    '''
    maxk = np.array([np.max([k1[i], k2[i], k3[i]]) for i in range(len(k1))])
    mink = np.array([np.min([k1[i], k2[i], k3[i]]) for i in range(len(k1))])
    if triangle == 'equilateral':       # only keep equilateral triangles
        triangle_index = np.where((k1 == k2) & (k2 == k3))
    elif triangle == 'acute':           # acute triangle
        triangle_index = np.where((k1**2 + k2**2 > k3**2) & (k2**2 + k3**2 > k1**2) & (k3**2 + k1**2 > k2**2))
    elif triangle == 'obtuse':          # obtuse triangle
        triangle_index = np.where((k1**2 + k2**2 < k3**2) | (k2**2 + k3**2 < k1**2) | (k3**2 + k1**2 < k2**2))
    elif triangle == 'extended':        # extended triangle
        triangle_index = np.where(maxk/mink > 3.0)
    
    return triangle_index




if __name__=='__main__': 
    #MeasureBisp(catalog='bigmd', correction='true', n_mocks=1, Ngrid=360, Nthreads=1)
    #notif(toaddr='changh20@gmail.com', subject='Bisp QPM true n_mock=100 Done') 
    #MeasureBisp(catalog='qpm', correction='upweight', n_mocks=range(51,101), Ngrid=360, Nthreads=10)
    #notif(toaddr='changh20@gmail.com', subject='Bisp QPM upweight n_mock=100 Done') 
    #Bisp_FcImpact(catalog='qpm')
    #Bisp_FcImpact(catalog='qpm', GM_klim=True)
    #Bisp_FcImpact(catalog='nseries', GM_klim=True)
    BispTriangle(catalog='nseries', GM_klim=True)
    #Bisp_CMASScomp(catalog='nseries', GM_klim=True, quadrupole=False)
    #Bisp_CMASScomp(catalog='qpm', GM_klim=False, quadrupole=False)
    #Bisp_CMASScomp(catalog='qpm', GM_klim=True, quadrupole=False)
