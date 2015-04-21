import numpy as np
import pylab as py
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import interpolate
import sys
import re


class powerspec: 
    def __init__(self): 
        self.catalog = ''
        self.label = ''
    def readFile(self, catalog='LasDamasGeo', catalog_param=[1, 'a'], corr='true', corr_param=[5.3, 1.0], \
            shotnoise='Irand'): 
        '''
        Read in P(k) for specified catalog and file 
        '''
        if catalog == 'LasDamasGeo': 
            if len(catalog_param) != 2:         # Specify naming convention 
                raise NameError('LasDamasGeo needs [number, "letter"]')
            if type(catalog_param[0]) is not int:
                raise NameError('LasDamasGeo needs [number, "letter"]')
            if type(catalog_param[1]) is not str:
                raise NameError('LasDamasGeo needs [number, "letter"]')
            file_dir = '/mount/riachuelo1/hahn/power/LasDamas/Geo/'
            file_prefix = 'power_sdssmock_gamma_lrgFull_zm_oriana'
            file_midfix = '_no.rdcz.'
            file_suffix = '.grid360.P020000.box3600'
            if corr == 'true':    # Specify fiber collision correction method
                file_corr = 'dat'
                corr_label = 'True' 
            elif (corr == 'delta') or (corr == 'upweight'): 
                file_corr = 'fibcoll.dat.delta'
                corr_label = 'Delta'
            elif corr == 'peak':  # Bestfit sigma_peak = 5.3
                if len(corr_param) != 2:
                    raise NameError('Peak correction: specify [sigma, fpeak]')
                if type(corr_param[0]) is not float:
                    raise NameError('Peak correction: specify [sigma, fpeak]')
                if type(corr_param[1]) is not float:
                    raise NameError('Peak correction: specify [sigma, fpeak]')
                file_corr = ''.join(['fibcoll.dat.peak.sigma', str(corr_param[0]), '.fpeak', str(corr_param[1])])
                corr_label = 'Peak(\sigma='+str(corr_param[0])+',f_{peak}='+str(corr_param[1])+')'
            if shotnoise not in ['Irand', 'Igal', 'Iflorian']: 
                raise NameError('Shot noise can only be [Irand, Igal, Iflorian]')
            elif shotnoise == 'Irand': 
                file_shotnoise = ''
            else: 
                file_shotnoise = '.'+shotnoise
                corr_label += ','+shotnoise 
            file_name = ''.join([file_dir, file_prefix, str(catalog_param[0]+100)[1:3], catalog_param[1], \
                    file_midfix, file_corr, file_shotnoise, file_suffix])    # get powerspectrum file name
        if catalog == 'LasDamasBox': 
            if len(catalog_param) != 1:
                raise NameError('LasDamasBox needs [number]')
            if type(catalog_param[0]) is not int:
                raise NameError('LasDamasBox needs [number]')
            file_dir = '/mount/riachuelo1/hahn/power/LasDamas/Box/'
            file_prefix = 'power_mock_gamma_lrg21p2_zmo_Oriana_1' 
            file_midfix = '_z0p342_fof_b0p2.zdist.'
            file_suffix = 'dat.grid360.P0.box3600'
            if corr=='true':    # Specify fiber collision correction method
                file_corr = ''
                corr_label = 'True'
            elif corr=='delta': 
                file_corr = 'fibcoll.delta.'
                corr_label = 'Delta'
            elif corr=='peak':  # Bestfit sigma_peak = 5.3
                if len(corr_param) != 1:
                    raise NameError('Peak correction: specify [sigma]')
                if type(corr_param[0]) is not float:
                    raise NameError('Peak correction: specify [sigma]')
                file_corr = ''.join(['fibcoll.peak.sigma.', str(corr_param[0]), '.'])
                corr_label = ''.join(['Peak(\sigma = ', str(corr_param[0]), ')'])
            file_name = ''.join([
                file_dir, file_prefix, str(catalog_param[0]+1000)[1:4], file_midfix, file_corr, file_suffix
                ])
        if catalog.lower() == 'pthalo':         # PTHALO
            if type(catalog_param) is not int:
                raise NameError('PTHalo needs [number]')
            file_dir = '/mount/riachuelo1/hahn/power/PTHalo/v11p0/'
            file_prefix = 'power_cmass_dr11_north_ir4'
            file_midfix = '.v11.0.wghtv.txt.'
            file_suffix = '.100randoms.grid960.P020000.box3600'
            if corr == 'true':    # Specify fiber collision correction method
                file_corr = 'noweight'
                corr_label = 'True' 
            elif (corr == 'delta') or (corr == 'upweight'): 
                file_corr = 'upweight'
                corr_label = 'Upweight'
            elif corr == 'peak':  # Bestfit sigma_peak = 5.3
                if len(corr_param) != 2:
                    raise NameError('Peak correction: specify [sigma, fpeak]')
                if type(corr_param[0]) is not float:
                    raise NameError('Peak correction: specify [sigma, fpeak]')
                if type(corr_param[1]) is not float:
                    raise NameError('Peak correction: specify [sigma, fpeak]')
                file_corr = ''.join([
                    'peak.sigma', str(corr_param[0]), '.fpeak', str(corr_param[1])])
                corr_label = 'Peak(\sigma='+str(corr_param[0])+',f_{peak}='+str(corr_param[1])+')'
            if shotnoise not in ['Irand', 'Igal', 'Iflorian', 'Icole', 'Ifkp']: 
                raise NameError('Shot noise can only be [Irand, Igal, Iflorian]')
            elif shotnoise == 'Irand': 
                file_shotnoise = ''
            elif shotnoise == 'Iflorian':  
                file_shotnoise = 'florian_shotnoise.'
                corr_label += ','+shotnoise 
            elif shotnoise == 'Icole': 
                file_shotnoise = 'cole_shotnoise.'
                corr_label += ','+shotnoise 
            elif shotnoise == 'Ifkp': 
                file_shotnoise = 'fkp_shotnoise.'
                corr_label += ','+shotnoise 
            file_name = ''.join([file_dir, file_prefix, str(catalog_param+1000)[1:4], 
                    file_midfix, file_shotnoise, file_corr, file_suffix])    # get powerspectrum file name
        if catalog == 'TilingMock':         # for tiling mock
            file_dir = '/mount/riachuelo1/hahn/power/tiling_mocks/'     # tiling mock P(k) directory 
            # file names are hardcoded since there is only one tiling mock  
            if corr == 'true': 
                file_name = ''.join([file_dir, 
                    'power_cmass-boss5003sector-icoll012.dat.grid360.P020000.box4000']) 
                corr_label = 'True'
            elif corr == 'delta': 
                file_name = ''.join([file_dir, 
                    'power_cmass-boss5003sector-icoll012.fibcoll.dat.delta.grid360.P020000.box4000'])
                corr_label = 'Delta'
            elif corr == 'peak':            # peak correction sigma hard-coded for now 
                file_name = ''.join([file_dir, 
                    'power_cmass-boss5003sector-icoll012.fibcoll.dat.peak.sigma5.98.fpeak1.0.grid360.P020000.box4000']) 
                corr_label = ''.join(['Peak(\sigma = 5.98)'])

        self.catalog = catalog
        self.label = ''.join([
            r"$P(k)_{\rm{", catalog, ";", corr_label, "}}$"
            ])
        self.filespec = catalog_param
        self.file = file_name   # store filename 
        data = np.loadtxt(self.file)    # read in file 
        self.k = data[:,0]
        self.Pk = data[:,1]

def average_powerspec(catalog='LasDamasGeo', catalog_paramrange=[range(1,41), ['a','b','c','d']], corr='true', corr_param=[5.3, 1.0], 
        shotnoise='Irand'): 
    '''
    Calculates the average P(k) for specified catalog and file specifier range
    '''
    count = 0 
    if catalog=='LasDamasGeo': 
        for first in catalog_paramrange[0]: 
            for second in catalog_paramrange[1]: 
                power = powerspec()
                power.readFile(catalog=catalog, catalog_param=[first, second], corr=corr, corr_param=corr_param, shotnoise=shotnoise)
                if (first==catalog_paramrange[0][0] and second==catalog_paramrange[1][0]): 
                    Pksum = np.array(power.Pk)
                    kvalues = np.array(power.k)
                else: 
                    Pksum += np.array(power.Pk)
                    if np.array_equal(kvalues, np.array(power.k))=='True': 
                        raise TypeError('k values do not match') 
                count += 1
    elif catalog=='LasDamasBox': 
        for first in catalog_paramrange[0]: 
            power = powerspec()
            power.readFile(catalog=catalog, catalog_param=[first], corr=corr, corr_param=corr_param)
            if (first==catalog_paramrange[0][0]): 
                Pksum = np.array(power.Pk)
                kvalues = np.array(power.k)
            else: 
                Pksum += np.array(power.Pk)
                if np.array_equal(kvalues, np.array(power.k))=='True': 
                    raise TypeError('k values do not match') 
            count += 1
    elif catalog == 'TilingMock': 
        power = powerspec()
        power.readFile(catalog=catalog, corr=corr, corr_param=corr_param) 
        Pksum = np.array(power.Pk)
        kvalues = np.array(power.k)
        count=1
    elif catalog.lower() == 'pthalo':               # for PTHalo mock catalogs
        for first in catalog_paramrange: 
            power = powerspec()
            power.readFile(catalog=catalog, catalog_param=first, corr=corr, corr_param=corr_param, shotnoise=shotnoise)
            if first == catalog_paramrange[0]: 
                Pksum = np.array(power.Pk)
                kvalues = np.array(power.k)
            else: 
                Pksum += np.array(power.Pk)
                if np.array_equal(kvalues, np.array(power.k))=='True': 
                    raise TypeError('k values do not match') 
            count += 1

    avg_label = ((power.label).split(r'\rm{'))[-1].rstrip('$').rstrip('}')
    avgPower = powerspec() 
    avgPower.catalog = catalog 
    avgPower.label = r"$\overline{P(k)}_{\rm{"+avg_label+"}}$"
    avgPower.filespec = catalog_paramrange
    avgPower.Pk = Pksum/float(count)
    avgPower.k  = kvalues
    return avgPower

def residual(arr1, arr2):
    if len(arr1) != len(arr2):
        raise TypeError("Input array lengths do not match.")
    else:
        resid = np.array([arr1[i]/arr2[i] for i in range(len(arr1))])
    return resid
      
# Plotting -----------------
def plot_powerspec_fibcolcorr_comparison(catalog='LasDamasGeo', catalog_paramrange=[range(1,41), ['a','b','c','d']],
        correction_method=['true', 'delta', 'peak'], corr_params=[[], [], [6.3, 1.0]], shotnoise='Irand', resid='False'): 
    '''
    Comparison of P(k)avg for different fibercollision correciton methods
    Make sure 'true' is first in the correction method list  
    '''
    prettyplot()

    corr_color = ['black', 'blue', 'red'] 
    fig_dir = '/home/users/hahn/research/figures/boss/fiber_collision/'
    for i_corr, correction in enumerate(correction_method): 
        # compute average[P(k)] for each correction method
        avg_Pk = average_powerspec(catalog=catalog, catalog_paramrange=catalog_paramrange, corr=correction, 
                corr_param=corr_params[i_corr], shotnoise=shotnoise)

        fig = plt.figure(1, figsize=(7, 8))
        sub = fig.add_subplot(111)

        if resid=='False':  # P(k) Comparison 
            if correction=='true': 
                sub.plot(avg_Pk.k, avg_Pk.Pk,\
                        color=corr_color[i_corr], label=avg_Pk.label)
            else: 
                sub.scatter(avg_Pk.k, avg_Pk.Pk,\
                        color=corr_color[i_corr], label=avg_Pk.label)
        else:               # P(k) residual comparison
            if correction=='true': 
                avg_Pk_true = avg_Pk
            else: 
                resid_label = ''.join([(avg_Pk.label).rstrip('$'), '/', (avg_Pk_true.label).lstrip('$')])
                sub.scatter(avg_Pk.k, residual(avg_Pk.Pk, avg_Pk_true.Pk),\
                        color=corr_color[i_corr], label=resid_label)

        if (correction=='peak'):        # specify correction label 
            if (catalog=='LasDamasGeo'): 
                corrlabel = ''.join([correction, '; \sigma = ', str(corr_params[i_corr][0]),\
                        r'f_{\rm{peak}}=', str(corr_params[i_corr][1])])
            elif (catalog=='LasDamasBox'): 
                corrlabel = ''.join([correction, '; \sigma = ', str(corr_params[i_corr][0])])
            elif (catalog=='TilingMock'): 
                corrlabel = ''.join([correction, '; \sigma = 5.98'])
            sigma_flag = '_'.join([str(corr_params[i_corr][i]) for i in range(len(corr_params[i_corr]))])
        else: 
            corrlabel = correction

    if shotnoise == 'Irand': 
        shotnoise_flag = ''
    elif shotnoise == 'Igal': 
        shotnoise_flag = '_Igal'

    if resid=='False': 
        if catalog=='LasDamasGeo': 
            ylimit = [10**2,10**5.5]
            ytext = 10**5.2
        elif catalog=='LasDamasBox': 
            ylimit = [10**1.5,10**4]
            ytext = 10**3.5
        elif catalog == 'TilingMock': 
            ylimit = [10**2,10**5.5]
            ytext = 10**5.2

        ylabel = 'P(k)'
        sub.set_yscale('log')
        fig_name = ''.join(['powerspec_', catalog, '_fibcoll_', '_'.join(correction_method), 
            sigma_flag, shotnoise_flag, '_comparison.png'])
    else: 
        ylimit = [0.9,1.2]
        ytext = 1.15
        ylabel = r"$\overline{P(k)}/\overline{P(k)}_{\rm{True}}$"
        fig_name = ''.join(['powerspec_', catalog, '_fibcoll_', '_'.join(correction_method), 
            sigma_flag, shotnoise_flag, '_residual_comparison.png'])

    if catalog=='LasDamasGeo': 
        sub.text(2.0*10**-3.0, ytext, str(np.max(catalog_paramrange[0]))+' Las Damas Geo Mocks', fontsize=15)

    sub.set_xscale('log')
    sub.set_xlim([10**-3,10**0])
    sub.set_ylim(ylimit)
    sub.set_xlabel('k', fontsize=20)
    sub.set_ylabel(ylabel, fontsize=20)
    sub.legend(loc='lower left',scatterpoints=1,prop={'size':14})
    sub.grid(True)
    fig.savefig(''.join([fig_dir, fig_name]))
    fig.clear()

def plot_powerspec_comparison(powerlist, resid='False'):
    '''
    Compares P(k) in catalog_param list of powerspec class objects.
    Input denominator first
    '''
    prettyplot() 
    colorlist = ['black', 'blue', 'red', 'green', 'yellow']  # Running out of color ideas...
    colors = colorlist[0:len(powerlist)] 
    file_label = ''
    for i_Pk, Pk in enumerate(powerlist): 
        fig = plt.figure(1, figsize=(7, 8))
        sub = fig.add_subplot(111)
        if resid == 'False':            # P(k) Comparison 
            sub.scatter(Pk.k, Pk.Pk,\
                    color=colors[i_Pk], label=Pk.label)
        else:                           # P(k) residual comparison
            if i_Pk == 0: 
                Pk1 = Pk  
            else: 
                resid_label = ''.join([(Pk.label).rstrip('$'), '/', (Pk1.label).lstrip('$')])
                sub.scatter(Pk.k, residual(Pk.Pk, Pk1.Pk),\
                        color=colors[i_Pk], label=resid_label)
        if (Pk.label).find('overline') != -1: 
            avg_flag = 'avg'
        else: 
            avg_flag = ''
        Pk_label = (((Pk.label).split(r'\rm{'))[-1].rstrip('$')).rstrip('}')
        Pk_label = re.sub(r"\(.+?\)", '', Pk_label) 
        Pk_label = ''.join([avg_flag, Pk_label.split(';')[0]]+(Pk_label.split(';')[1]).split(','))
        file_label += '_'+Pk_label 
    if resid=='False': 
        ylimit = [10**2,10**5.5]
        ylabel = 'P(k)'
        sub.set_yscale('log')
        fig_name = ''.join(['powerspec_', file_label, '_comparison.png'])
    else: 
        ylimit = [0.9,1.1]
        ylabel = resid_label 
        fig_name = ''.join(['powerspec_', file_label, '_residual_comparison.png'])
    sub.set_xscale('log')
    sub.set_xlim([10**-3,10**0])
    sub.set_ylim(ylimit)
    sub.set_xlabel('k',fontsize=20)
    sub.set_ylabel(ylabel,fontsize=20)
    sub.legend(loc='lower left',prop={'size':14}, markerscale=2.5, scatterpoints=1)

    sub.grid(True)
    fig_dir = '/home/users/hahn/research/figures/boss/fiber_collision/'
    fig.savefig(''.join([fig_dir, fig_name]), bbox_inches=0)
    fig.clear()

# Analysis -----------------
def Pk_shotnoise_comparison(catalog='LasDamasGeo'): 
    '''
    Compares P(k) derived from different shotnoise computation methods
    '''
    if catalog == 'LasDamasGeo': 
        # For single P(k)
        pk_true = powerspec()              # import reference true P(k) with no Fiber Collisions
        pk_true.readFile(catalog=catalog, corr='true')
        pk_upweight = powerspec()           # import upweight with FKP shot noise correction 
        pk_upweight.readFile(catalog=catalog, corr='delta')
        pk_florian = powerspec()            # import upweight with Florian shot nosie correction 
        pk_florian.readFile(catalog=catalog, corr='delta', shotnoise='Iflorian')
        pk_gal = powerspec()            # import upweight with galaxy shot nosie correction 
        pk_gal.readFile(catalog=catalog, corr='delta', shotnoise='Igal')
        avgpk_true = average_powerspec(catalog=catalog, catalog_paramrange=[range(1,6), ['a','b','c','d']], corr='true')
        avgpk_upweight = average_powerspec(catalog=catalog, catalog_paramrange=[range(1,6), ['a','b','c','d']], corr='upweight')
        avgpk_gal = average_powerspec(catalog=catalog, catalog_paramrange=[range(1,6), ['a','b','c','d']], 
                corr='upweight', shotnoise='Igal')
        avgpk_florian = average_powerspec(catalog=catalog, catalog_paramrange=[range(1,6), ['a','b','c','d']], 
                corr='upweight', shotnoise='Iflorian')

        plot_powerspec_comparison([pk_true, pk_upweight, pk_gal, pk_florian], resid='False')
        plot_powerspec_comparison([pk_true, pk_upweight, pk_gal, pk_florian], resid='True')

        plot_powerspec_comparison([avgpk_true, avgpk_upweight, avgpk_gal, avgpk_florian], resid='False')
        plot_powerspec_comparison([avgpk_true, avgpk_upweight, avgpk_gal, avgpk_florian], resid='True')
    elif catalog.lower() == 'pthalo': 
        # import all different shotnoise method P(k)
        pk_true = powerspec()                           # import reference true P(k) with no Fiber Collisions
        pk_true.readFile(catalog=catalog, catalog_param=1, corr='true', shotnoise='Ifkp')
        pk_upweight = powerspec()                       # import upweight with ngalsys FKP shot noise correction 
        pk_upweight.readFile(catalog=catalog, catalog_param=1, corr='delta')
        pk_florian = powerspec()                        # import upweight with Florian shot nosie correction 
        pk_florian.readFile(catalog=catalog, catalog_param=1, corr='delta', shotnoise='Iflorian')
        pk_cole = powerspec()                           # import upweight with Cole shot nosie correction 
        pk_cole.readFile(catalog=catalog, catalog_param=1, corr='delta', shotnoise='Icole')
        pk_fkp = powerspec()                            # import upweight with original FKP shot nosie correction 
        pk_fkp.readFile(catalog=catalog, catalog_param=1, corr='delta', shotnoise='Ifkp')

        plot_powerspec_comparison([pk_true, pk_upweight, pk_florian, pk_cole, pk_fkp], resid='False')
        plot_powerspec_comparison([pk_true, pk_upweight, pk_florian, pk_cole, pk_fkp], resid='True')

        #avgpk_true = average_powerspec(catalog=catalog, catalog_paramrange=range(1,11), corr='true')
        #avgpk_upweight = average_powerspec(catalog=catalog, catalog_paramrange=range(1,11), corr='upweight')
        #avgpk_florian = average_powerspec(catalog=catalog, catalog_paramrange=range(1,11), 
        #        corr='upweight', shotnoise='Iflorian')

        #plot_powerspec_comparison([avgpk_true, avgpk_upweight, avgpk_florian], resid='False')
        #plot_powerspec_comparison([avgpk_true, avgpk_upweight, avgpk_florian], resid='True')
        #
        #plot_powerspec_comparison([avgpk_upweight, avgpk_florian], resid='True')

def LasDamasGeo_fibcoll_comparison():
    '''
    Compare average P(k) for LasDamasGeo Mock Catalogs corrected using True, Delta, and Peak fiber collision correction methods
    '''
    #plot_powerspec_fibcolcorr_comparison(catalog='LasDamasGeo', catalog_paramrange=[range(1, 41), ['a', 'b', 'c', 'd']],
    #        correction_method=['true', 'delta', 'peak'], corr_params=[[], [], [5.3, 1.0]])
    #plot_powerspec_fibcolcorr_comparison(catalog='LasDamasGeo', catalog_paramrange=[range(1, 41), ['a', 'b', 'c', 'd']],
    #        correction_method=['true', 'delta', 'peak'], corr_params=[[], [], [5.3, 1.0]], resid='True')
    plot_powerspec_fibcolcorr_comparison(catalog='LasDamasGeo', catalog_paramrange=[range(1, 41), ['a', 'b', 'c', 'd']],
            correction_method=['true', 'delta', 'peak'], corr_params=[[], [], [5.3, 1.0]], shotnoise='Igal')
    plot_powerspec_fibcolcorr_comparison(catalog='LasDamasGeo', catalog_paramrange=[range(1, 41), ['a', 'b', 'c', 'd']],
            correction_method=['true', 'delta', 'peak'], corr_params=[[], [], [5.3, 1.0]], shotnoise='Igal', resid='True')

if __name__=='__main__': 
    #Pk_shotnoise_comparison(catalog='pthalo')
    #LasDamasGeo_fibcoll_comparison()
    plot_powerspec_fibcolcorr_comparison(catalog='LasDamasGeo', catalog_paramrange=[range(1,2), ['a','b','c','d']], 
            correction_method=['true', 'peak'], corr_params=[[], [5.3, 0.1]])
    plot_powerspec_fibcolcorr_comparison(catalog='LasDamasGeo', catalog_paramrange=[range(1,2), ['a','b','c','d']], 
            correction_method=['true', 'peak'], corr_params=[[], [5.3, 0.1]], resid='True')
    #plot_powerspec_fibcolcorr_comparison(catalog='TilingMock', \
    #        correction_method=['true', 'delta', 'peak'], corr_params=[[], [], [5.98]], resid='True')
    #plot_powerspec_fibcolcorr_comparison(catalog='PTHalo', \
    #        correction_method=['true', 'delta', 'peak'], corr_params=[[], [], [5.98]])
    #plot_powerspec_fibcolcorr_comparison(catalog='TilingMock', \
    #        correction_method=['true', 'delta', 'peak'], corr_params=[[], [], [5.98]], resid='True')
    #plot_powerspec_fibcolcorr_comparison(catalog='LasDamasBox', catalog_paramrange=[range(1,2)], 
    #        correction_method=['true', 'delta', 'peak'], corr_params=[[], [], [5.0]], resid='True')
    #plot_powerspec_fibcolcorr_comparison(catalog='LasDamasBox', catalog_paramrange=[range(1,2)], \
    #        correction_method=['true', 'delta', 'peak'], corr_params=[[], [], [7.6]])
    #plot_powerspec_fibcolcorr_comparison(catalog='LasDamasBox', catalog_paramrange=[range(1,2)], 
    #        correction_method=['true', 'delta', 'peak'], corr_params=[[], [], [7.6]], resid='True')
