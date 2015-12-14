'''

tests of Line of Sight Displacement Codes

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
    ''' compute dlos using cosmolopy and compare dLOS calculated from IDL  
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

