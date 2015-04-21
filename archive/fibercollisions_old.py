import numpy as np
import pylab as py
import math as m
import matplotlib.pyplot as plt
from matplotlib import rc
import os.path
import subprocess
import cosmolopy as cosmos

# Classes ------------------------------------------------------------
class galaxy_data: 
    def __init__(self, DorR='data', catalog='LasDamasGeo', catalog_param={'n_mock': [1,'a']}, 
            corr='peak', corr_param={'sigma':5.3, 'fpeak':0.5}, readdata=True): 
        '''
        Read in the data or random file and store all the appropriate values  
        '''
        if catalog.lower() == 'lasdamasgeo': 
            # For data 
            if DorR == 'data': 
                # columns that this catalog data will have  
                catalog_columns = ['ra', 'dec', 'z', 'weight']
                self.columns = catalog_columns

                # True LasDamasGeo (No FiberCollisions) 
                if corr.lower() == 'true': 
                    file_name = '/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/sdssmock_gamma_lrgFull_zm_oriana'+\
                            str("%02d" % catalog_param['n_mock'][0])+catalog_param['n_mock'][1]+'_no.rdcz.dat'

                    if readdata == True: 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])         # ra, dec, CZ

                        # assign to class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            if catalog_column == 'z': 
                                column_data = file_data[i_col]/299800.0     # divide by speed of light
                            elif catalog_column == 'weight':                    
                                column_data = np.array([1.0 for i_data in range(len(file_data[:,0]))]) # no weights for true (all 1) 
                            else: 
                                column_data = file_data[i_col]
                            # assign to class
                            setattr(self, catalog_column, column_data) 

                # Fibercollisions Corrected by Upweight correction 
                elif corr.lower() == 'upweight': 
                    # correction parameters not needed
                    file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+\
                            str("%02d" % catalog_param['n_mock'][0])+catalog_param['n_mock'][1]+'_no.rdcz.fibcoll.dat'

                    if readdata == True: 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights 

                        # assign to class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 

                # Fibercollisions Corrected by peak correction 
                elif (corr.lower()  == 'peak') or (corr.lower() == 'peaknbar'): 
                    corr_str = '.peak.sigma'+str(corr_param['sigma'])+'.fpeak'+str(corr_param['fpeak'])

                    file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+\
                            str("%02d" % catalog_param['n_mock'][0])+catalog_param['n_mock'][1]+'_no.rdcz.fibcoll.dat'+corr_str
                    
                    if readdata == True: 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights 

                        # assign to class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 
                
                elif corr.lower()  == 'peaktest': 
                    corr_str = '.peaktest.sigma'+str(corr_param['sigma'])+'.fpeak'+str(corr_param['fpeak'])

                    file_name = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull_zm_oriana'+\
                            str("%02d" % catalog_param['n_mock'][0])+catalog_param['n_mock'][1]+'_no.rdcz.fibcoll.dat'+corr_str
                    
                    if readdata == True: 
                        file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2,3])         # ra, dec, z, weights 

                        # assign to class
                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            setattr(self, catalog_column, column_data) 

            # For Random 
            else: 
                # columns that this catalog random will have
                catalog_columns = ['ra', 'dec', 'z'] 
                self.columns = catalog_columns
                
                # random for True LasDamasGeo data 
                if corr.lower() == 'true': 
                    file_name = '/mount/chichipio2/rs123/MOCKS/randoms/sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat'
                    file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2])

                    for i_col, catalog_column in enumerate(catalog_columns): 
                        if catalog_column == 'z': 
                            column_data = file_data[i_col]/299800.0
                        else: 
                            column_data = file_data[i_col]
                        # assign to class
                        setattr(self, catalog_column, column_data)
                
                # random corrected for fiber collision upweight correction: 
                elif corr.lower() == 'upweight': 
                    corr_str = '.upweight'

                    file_name = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/', 
                        'sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat.allmocks', corr_str]) 
                    file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2]) 

                    for i_col, catalog_column in enumerate(catalog_columns): 
                        column_data = file_data[i_col]
                        # assign to class
                        setattr(self, catalog_column, column_data)

                # random corrected for fiber collision peak correction: 
                elif (corr.lower() == 'peak') or (corr.lower() == 'peaknbar'): 
                    # currently hardcoded a bit
                    corr_str = '.peak.sigma'+str(corr_param['sigma'])+'.fpeak'+str(corr_param['fpeak'])

                    file_name = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/',
                        'sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat.allmocks', corr_str])

                    if readdata == True:  
                        # check if file exists 
                        if os.path.isfile(file_name) == True: 
                            file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2]) 
                        else:       # if it doesn't generate the corrected randoms
                            print 'BUILDING CORRECTED RANDOM FILE ', file_name
                            build_corrected_randoms(catalog=catalog, 
                                    catalog_param_list={'n_mock':[[i,j] for i in range(1,41) for j in ['a', 'b', 'c', 'd']]}, 
                                    corr=corr.lower(), corr_param={'sigma':5.3, 'fpeak':corr_param['fpeak']}, sanitycheck=False)
                            file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2]) 

                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            # assign to class
                            setattr(self, catalog_column, column_data)
                
                # random corrected for fiber collision peak test correction: 
                elif corr.lower() == 'peaktest': 
                    # currently hardcoded a bit
                    corr_str = '.peaktest.sigma'+str(corr_param['sigma'])+'.fpeak'+str(corr_param['fpeak'])

                    file_name = ''.join(['/mount/riachuelo1/hahn/data/LasDamas/Geo/',
                        'sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat.allmocks', corr_str])

                    if readdata == True:  
                        # check if file exists 
                        if os.path.isfile(file_name) == True: 
                            file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2]) 
                        else:       # if it doesn't generate the corrected randoms
                            print 'BUILDING CORRECTED RANDOM FILE ', file_name
                            build_corrected_randoms(catalog=catalog, 
                                    catalog_param_list={'n_mock':[[i,j] for i in range(1,41) for j in ['a', 'b', 'c', 'd']]}, 
                                    corr=corr.lower(), corr_param={'sigma':5.3, 'fpeak':corr_param['fpeak']}, sanitycheck=False)
                            file_data = np.loadtxt(file_name, unpack=True, usecols=[0,1,2]) 

                        for i_col, catalog_column in enumerate(catalog_columns): 
                            column_data = file_data[i_col]
                            # assign to class
                            setattr(self, catalog_column, column_data)

            # save filename (for accounting) 
            self.file_name = file_name 

class spec: 
    def __init__(self, spectrum='power', catalog='tilingmock', filespec={'version':'v11p0', 'n_mock': 1, 'Nrandom':100}, 
            fibcolcorr='true', corrspec={'sigma':5.3,'fpeak':1.0}): 
        '''
        bispectrum file 
        specify catalog, version, mock file number, file specifications (e.g. Nrandom), 
        fiber collision correction method, correction specifications (e.g. sigma, fpeak)
        '''
        # powerspectrum or bispectrum
        self.spectrum = spectrum
        if spectrum == 'power':                         
            spec_dir_flag = 'power'
            spec_file_flag = 'power_'
        elif spectrum == 'bispec': 
            spec_dir_flag = 'bispec'
            spec_file_flag = 'bisp_'

        # Tiling Mock 
        if catalog.lower() == 'tilingmock':  
            file_dir = '/mount/riachuelo1/hahn/'+spec_dir_flag+'/tiling_mocks/'      
            file_prefix = spec_file_flag+'cmass-boss5003sector-icoll012.'
            if spectrum == 'bispec': 
                file_suffix = '.grid360.nmax.nstep3.P020000.box4000'
            elif spectrum == 'power': 
                file_suffix = '.grid360.P020000.box4000'

            if fibcolcorr == 'peak': 
                # e.g. power_cmass-boss5003sector-icoll012.fibcoll.dat.peak.sigma5.98.fpeak1.0.grid360.P020000.box4000
                file_corr = 'fibcoll.dat.peak.sigma5.98.fpeak1.0'
            elif fibcolcorr == 'delta': 
                # e.g. power_cmass-boss5003sector-icoll012.fibcoll.dat.delta.grid360.P020000.box4000
                file_corr = 'fibcoll.dat.delta'
            else: 
                file_corr = 'dat'
            self.scale = 4000

        # PT Halo Mock 
        elif catalog.lower() == 'pthalo':
            file_dir = ''.join(['/mount/riachuelo1/hahn/', spec_dir_flag,'/PTHalo/', filespec['version'].lower(), '/']) 
            
            if filespec['version'].lower() == 'v11p0':    # specify version 
                # e.g. bisp_cmass_dr11_north_ir4001.v11.0.wghtv.txt.upweight.100randoms.grid360.nmax40.nstep3.P020000.box3600
                # e.g. power_cmass_dr11_north_ir4010.v11.0.wghtv.txt.upweight.100randoms.grid360.P020000.box3600
                file_prefix = ''.join([spec_file_flag, 'cmass_dr11_north_ir4', str(1000+filespec['n_mock'])[1:4]])  
                if spectrum == 'power': 
                    file_suffix = ''.join(['.', str(filespec['Nrandom']), 'randoms.grid360.P020000.box3600']) 
                elif spectrum == 'bispec': 
                    file_suffix = ''.join(['.', str(filespec['Nrandom']), 'randoms.grid360.nmax40.nstep3.P020000.box3600']) 
                
                if fibcolcorr == 'peak':        # specify fiber collision method 
                    file_corr = ''.join(['.v11.0.wghtv.txt.peak.sigma', str(corrspec['sigma']), '.fpeak', str(corrspec['fpeak'])])
                elif fibcolcorr == 'delta': 
                    file_corr = '.v11.0.wghtv.txt.upweight'
                elif fibcolcorr == 'wboss': 
                    file_corr = '.v11.0.wghtv.txt.wboss'
                else: 
                    file_corr = '.v11.0.wghtv.txt.noweights'
            self.scale = 3600   # Mpc (Survey scale) 

        # LasDamas Geomtry Mocks
        elif catalog.lower() == "lasdamasgeo": 
            # Las Damas Directory
            file_dir = ''.join(['/mount/riachuelo1/hahn/', spec_dir_flag, '/LasDamas/Geo/'])

            # e.g. power_sdssmock_gamma_lrgFull_zm_oriana19a_no.rdcz.dat.grid360.P020000.box3600
            # file beginning
            file_prefix = spec_file_flag+'sdssmock_gamma_lrgFull_zm_oriana'+str(100+filespec['n_mock'][0])[1:3]+filespec['n_mock'][1]+\
                    '_no.rdcz.'

            # file ending  
            if spectrum == 'bispec': 
                file_suffix = '.grid360.nmax.nstep3.P020000.box3600'
            elif spectrum == 'power': 
                file_suffix = '.grid360.P020000.box3600'
            
            # correction specifier 
            if fibcolcorr == 'peak': 
                # e.g. power_sdssmock_gamma_lrgFull_zm_oriana01c_no.rdcz.fibcoll.dat.peak.sigma5.3.fpeak0.7.grid360.P020000.box3600
                file_corr = 'fibcoll.dat.peak.sigma'+str(corrspec['sigma'])+'.fpeak'+str(corrspec['fpeak'])+'.corrnbar'
            elif fibcolcorr == 'delta': 
                # e.g. power_cmass-boss5003sector-icoll012.fibcoll.dat.delta.grid360.P020000.box4000
                file_corr = 'fibcoll.dat.delta'
            elif fibcolcorr == 'randdisptest': 
                file_corr = 'randdisptest.dat.randdisttest'
            else: 
                file_corr = 'dat'

            # survey scale 
            self.scale = 3600 
        
        # combine to make file
        self.file = ''.join([file_dir, file_prefix, file_corr, file_suffix]) # combine file parts 
        print self.file

    def readfile(self): 
        '''
        Read power/bi-spectrum file and import values of interest  
        '''
        file = np.loadtxt(self.file) 
        if self.spectrum == 'power':
            # columns of data to be imported
            self.columns = ['k', 'Pk']
            self.k  = file[:,0]
            self.Pk = file[:,1]
        elif self.spectrum == 'bispec': 
            # columns of data to be imported
            self.columns = ['kfund', 'i_triangles', 'k1', 'k2', 'k3', 'Pk1', 'Pk2', 'Pk3', 
                    'Bk', 'Q', 'avgk', 'kmax']

            k_fund = (2.0*m.pi)/np.float(self.scale)
            self.kfund = k_fund     # fundamental k 
            # import values from file  
            self.i_triangle = range(len(file[:,0]))
            self.k1 = k_fund*file[:,0]
            self.k2 = k_fund*file[:,1]
            self.k3 = k_fund*file[:,2]
            self.Pk1 = file[:,3]
            self.Pk2 = file[:,4]
            self.Pk3 = file[:,5]
            self.Bk = file[:,6]
            self.Q = file[:,7]
            self.avgk = np.array([np.mean([self.k1[i], self.k2[i], self.k3[i]]) for i in range(len(self.k1))])
            self.kmax = np.array([np.max([self.k1[i], self.k2[i], self.k3[i]]) for i in range(len(self.k1))])

class nbar: 
    def __init__(self, catalog='lasdamasgeo', corr='peaknbar', corr_param={'sigma':5.3, 'fpeak':0.1}): 
        '''
        read in appropriate nbar(z) file
        '''
        # set zlow, zmid, zhigh (consistent with CMASS nbar(z) files) 
        self.zlow = np.array([0.005*np.float(i) for i in range(201)])
        self.zhigh = self.zlow+0.005 
        self.zmid = self.zlow+0.0025 

        # store catalog info 
        self.catalog = catalog 
        self.corr = corr
        self.corr_param = corr_param

        corr_param_str = ''
        for key in self.corr_param.keys(): 
            corr_param_str = corr_param_str+key+str(self.corr_param[key])
        
        if self.catalog.lower() == 'lasdamasgeo': 
            nbar_file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'

        self.file_name = nbar_file_dir+'-'.join(['nbar', self.catalog.lower(), self.corr, corr_param_str])+'.dat'

        # for LasDamasGeo
        if catalog.lower() == 'lasdamasgeo': 
            ldg_nbar = 9.44451*10.0**-5         # constant nbar(z) value for true

            # for no correction, constant nbar(z)=9.44451e-5 
            if corr.lower() == 'true': 
                self.nbar = np.array([ldg_nbar for i in range(len(self.zlow))]) 

            # for upweight correction
            elif corr.lower() == 'upweight': 
                if os.path.isfile(self.file_name) == True: 
                    self.nbar = np.loadtxt(self.file_name, unpack=True, usecols=[3])
                else: 
                    # check that the ngal files exist
                    corr_ngal_file = get_nbar_ngal_file(DorR='random', catalog=catalog, corr='upweight')
                    if os.path.isfile(corr_ngal_file) == False: 
                        print corr_ngal_file, ' DOES NOT EXIST YET!'
                        write_nbar_ngal(DorR='random', catalog=catalog, corr='upweight') 

                    # read corrected nbar_ngal
                    corr_ngal = np.loadtxt(corr_ngal_file, unpack=True, usecols=[3]) 
       
                    # check if it exists
                    true_ngal_file = get_nbar_ngal_file(DorR='random', catalog=catalog, corr='true') 
                    if os.path.isfile(true_ngal_file) == False: 
                        print true_ngal_file, ' DOES NOT EXIST YET!'
                        write_nbar_ngal(DorR='random', catalog=catalog, corr='true', corr_param=corr_param)
                    
                    # read true nbar_ngal
                    true_ngal = np.loadtxt(true_ngal_file, unpack=True, usecols=[3])

                    # use corrected nbar_ngal file to determine corrected nbar(z) 
                    self.nbar = np.zeros(len(self.zlow))
                    for i in range(len(true_ngal)): 
                        if true_ngal[i] != 0: 
                            self.nbar[i] = ldg_nbar*(corr_ngal[i]/true_ngal[i])
                        else: 
                            self.nbar[i] = 0.0
                    self.writenbar()

            # for peak nbar correction, scale using ngal 
            elif (corr.lower() == 'peak') or (corr.lower() == 'peaknbar'): 
                if os.path.isfile(self.file_name) == True: 
                    self.nbar = np.loadtxt(self.file_name, unpack=True, usecols=[3])
                else: 
                    # check that the ngal files exist
                    corr_ngal_file = get_nbar_ngal_file(DorR='random', catalog=catalog, corr='peaknbar', corr_param=corr_param)
                    if os.path.isfile(corr_ngal_file) == False: 
                        print corr_ngal_file, ' DOES NOT EXIST YET!'
                        write_nbar_ngal(DorR='random', catalog=catalog, corr='peaknbar', corr_param=corr_param)

                    # read corrected nbar_ngal
                    corr_ngal = np.loadtxt(corr_ngal_file, unpack=True, usecols=[3]) 
       
                    # check if it exists
                    true_ngal_file = get_nbar_ngal_file(DorR='random', catalog=catalog, corr='true') 
                    if os.path.isfile(true_ngal_file) == False: 
                        print true_ngal_file, ' DOES NOT EXIST YET!'
                        write_nbar_ngal(DorR='random', catalog=catalog, corr='true', corr_param=corr_param)
                    
                    # read true nbar_ngal
                    true_ngal = np.loadtxt(true_ngal_file, unpack=True, usecols=[3])

                    # use corrected nbar_ngal file to determine corrected nbar(z) 
                    self.nbar = np.zeros(len(self.zlow))
                    for i in range(len(true_ngal)): 
                        if true_ngal[i] != 0: 
                            self.nbar[i] = ldg_nbar*(corr_ngal[i]/true_ngal[i])
                        else: 
                            self.nbar[i] = 0.0
                    self.writenbar()
            elif corr.lower() == 'peaktest': 
                if os.path.isfile(self.file_name) == True: 
                    self.nbar = np.loadtxt(self.file_name, unpack=True, usecols=[3])
                else: 
                    # check that the ngal files exist
                    corr_ngal_file = get_nbar_ngal_file(DorR='random', catalog=catalog, corr=corr, corr_param=corr_param)
                    if os.path.isfile(corr_ngal_file) == False: 
                        print corr_ngal_file, ' DOES NOT EXIST YET!'
                        write_nbar_ngal(DorR='random', catalog=catalog, corr='peaknbar', corr_param=corr_param)

                    # read corrected nbar_ngal
                    corr_ngal = np.loadtxt(corr_ngal_file, unpack=True, usecols=[3]) 
       
                    # check if it exists
                    true_ngal_file = get_nbar_ngal_file(DorR='random', catalog=catalog, corr='true') 
                    if os.path.isfile(true_ngal_file) == False: 
                        print true_ngal_file, ' DOES NOT EXIST YET!'
                        write_nbar_ngal(DorR='random', catalog=catalog, corr='true', corr_param=corr_param)
                    
                    # read true nbar_ngal
                    true_ngal = np.loadtxt(true_ngal_file, unpack=True, usecols=[3])

                    # use corrected nbar_ngal file to determine corrected nbar(z) 
                    self.nbar = np.zeros(len(self.zlow))
                    for i in range(len(true_ngal)): 
                        if true_ngal[i] != 0: 
                            self.nbar[i] = ldg_nbar*(corr_ngal[i]/true_ngal[i])
                        else: 
                            self.nbar[i] = 0.0
                    self.writenbar()
            else: 
                raise NameError("Not yet Coded!") 
        else: 
            raise NameError("Not yet Coded!") 

    def writenbar(self): 
        '''
        Write class object nbar to ASCII file 
        '''
        np.savetxt(self.file_name, 
                np.c_[self.zmid, self.zlow, self.zhigh, self.nbar], 
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5e'], delimiter='\t')

# Functions -----------------------------------------------------------------
def get_nbar_ngal_file(DorR='data', catalog='lasdamasgeo', catalog_param={'n_mock':[1,'a']},
        corr='peaknbar', corr_param={'sigma': 5.3, 'fpeak': 0.1}):
    '''
    get nbar_ngal file name
    '''
    file_prefix = 'nbar-ngal-'
    
    if catalog.lower() == 'lasdamasgeo': 
        catalog_str = catalog.lower()
        
        file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
        
        # specify data or random
        DorR_str = DorR.lower()

        # if data, specify specific catalog #  
        if DorR.lower() == 'data':
            catalog_str = catalog_str+'-'+str(catalog_param['n_mock'][0])+catalog_param['n_mock'][1]
    
        # specify correction
        if corr.lower() == 'true':
            corr_str = 'true'
        elif (corr.lower() == 'peak') or (corr.lower() == 'peaknbar'):
            corr_str = 'peak.sigma'+str(corr_param['sigma'])+'fpeak'+str(corr_param['fpeak'])
        elif corr.lower() == 'peaktest': 
            corr_str = 'peaktest.sigma'+str(corr_param['sigma'])+'fpeak'+str(corr_param['fpeak'])
        elif corr.lower() == 'upweight':
            corr_str = 'fibcoll.upweight'
    else: 
        raise NameError('not yet coded!')

    # combine to form filename  
    file_name = ''.join([file_dir, file_prefix, catalog_str, '-', DorR_str, '-', corr_str, '.dat'])
    return file_name

def write_nbar_ngal(DorR='data', catalog='lasdamasgeo', catalog_param={'n_mock': [1, 'a']}, 
        corr='peaknbar', corr_param={'sigma':5.3, 'fpeak':0.1}): 
    '''
    write ngal values for nbar(z) redshift bins to a file so it doesn't have to be repeated
    '''
    
    zcen, zlow, zhigh, shell_vol = np.loadtxt('/mount/riachuelo1/hahn/data/nbar-junk.dat', unpack=True, usecols=[0,1,2,5])
    z_values =[zcen, zlow, zhigh, shell_vol]

    if DorR.lower() == 'random':
        # True Random catalog
        if corr.lower() == 'true':
            # Hardcoded for LasDamasGeo
            if catalog.lower() == 'lasdamasgeo':
                true_rand_file = '/mount/chichipio2/rs123/MOCKS/randoms/sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat'
            else:
                raise NameError('Only LasDamasGeo coded so far')

            # read in files
            true_rand_ra, true_rand_dec, true_rand_cz = np.loadtxt(true_rand_file, unpack=True, usecols=[0,1,2])

            # get z values and the weights
            z_dist = true_rand_cz/299800.0
            z_weights = np.array([1 for i in range(len(z_dist))])
        # upweight corrected Random catalog
        elif corr.lower() == 'upweight':
            if catalog.lower() == 'lasdamasgeo':
                corr_rand = galaxy_data(DorR='random', catalog='LasDamasGeo', corr=corr, 
                        corr_param=corr_param) 
            else:
                raise NameError('Only LasDamasGeo coded so far')
            # read corrected random file
            corr_rand_ra = corr_rand.ra
            corr_rand_dec = corr_rand.dec
            corr_rand_z = corr_rand.z

            # get z values and weights
            z_dist = corr_rand_z
            z_weights = np.array([1 for i in range(len(z_dist))])

        # Peak corrected Random catalog
        elif (corr.lower() == 'peak') or (corr.lower() == 'peaknbar'):
            if catalog.lower() == 'lasdamasgeo':
                corr_rand = galaxy_data(DorR='random', catalog='LasDamasGeo', corr=corr, 
                        corr_param=corr_param) 
            else:
                raise NameError('Only LasDamasGeo coded so far')
            # read corrected random file
            corr_rand_ra = corr_rand.ra
            corr_rand_dec = corr_rand.dec
            corr_rand_z = corr_rand.z

            # get z values and weights
            z_dist = corr_rand_z
            z_weights = np.array([1 for i in range(len(z_dist))])
        # Peak Test corrected Random catalog
        elif corr.lower() == 'peaktest':
            if catalog.lower() == 'lasdamasgeo':
                corr_rand = galaxy_data(DorR='random', catalog='LasDamasGeo', corr=corr, 
                        corr_param=corr_param) 
            else:
                raise NameError('Only LasDamasGeo coded so far')
            # read corrected random file
            corr_rand_ra = corr_rand.ra
            corr_rand_dec = corr_rand.dec
            corr_rand_z = corr_rand.z

            # get z values and weights
            z_dist = corr_rand_z
            z_weights = np.array([1 for i in range(len(z_dist))])

    # data catalog
    elif DorR.lower() == 'data':
        # For true data
        if corr.lower() == 'true':
            if catalog.lower() == 'lasdamasgeo':
                true_file = ''.join(['/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/',
                    'sdssmock_gamma_lrgFull_zm_oriana', str("%02d" % catalog_param['n_mock'][0]), catalog_param['n_mock'][1],
                    '_no.rdcz.dat'])
                # read true file
                true_ra, true_dec, true_cz = np.loadtxt(true_file, unpack=True, usecols=[0,1,2])

                # get z values and weights 
                z_dist = true_cz/299800.0
                z_weights = np.array([1 for i in range(len(z_dist))])

        # For upweight corrected data
        elif corr.lower() == 'upweight':
            if catalog.lower() == 'lasdamasgeo':
                corr_data = galaxy_data(DorR='data', catalog=catalog, catalog_param=catalog_param, corr=corr, 
                        corr_param=corr_param )
            else: 
                raise NameError('Only LasDamasGeo coded so far')
            # read-in corrected file
            corr_ra = corr_data.ra 
            corr_dec = corr_data.dec
            corr_z = corr_data.z
            corr_w = corr_data.weight

            # get z values and weights
            z_dist = corr_z
            z_weights = corr_w

        # For peak corrected data
        elif (corr.lower() == 'peak') or (corr.lower() == 'peaknbar'):
            if catalog.lower() == 'lasdamasgeo':
                corr_data = galaxy_data(DorR='data', catalog=catalog, catalog_param=catalog_param, corr=corr, 
                        corr_param=corr_param )
            else: 
                raise NameError('Only LasDamasGeo coded so far')
            # read-in corrected file
            corr_ra = corr_data.ra 
            corr_dec = corr_data.dec
            corr_z = corr_data.z
            corr_w = corr_data.weight

            # get z values and weights
            z_dist = corr_z
            z_weights = corr_w
        # For peak test corrected data
        elif corr.lower() == 'peaktest':
            if catalog.lower() == 'lasdamasgeo':
                corr_data = galaxy_data(DorR='data', catalog=catalog, catalog_param=catalog_param, corr=corr, 
                        corr_param=corr_param )
            else: 
                raise NameError('Only LasDamasGeo coded so far')
            # read-in corrected file
            corr_ra = corr_data.ra 
            corr_dec = corr_data.dec
            corr_z = corr_data.z
            corr_w = corr_data.weight

            # get z values and weights
            z_dist = corr_z
            z_weights = corr_w
        else: 
            raise NameError('Not yet corrected!')

    nbar_ngal = np.zeros(len(z_values[0]))
    for i_z, zmid in enumerate(z_values[0]):
        zlim = (z_dist >= (z_values[1])[i_z]) & (z_dist < (z_values[2])[i_z])
        nbar_ngal[i_z] = np.sum(z_weights[zlim])
       
        # write nbar_ngal data to ask ascii file
        nbar_ngal_file = get_nbar_ngal_file(DorR=DorR, catalog=catalog, catalog_param=catalog_param,
                corr=corr, corr_param=corr_param)
        np.savetxt(nbar_ngal_file,            
                np.c_[z_values[0], z_values[1], z_values[2], nbar_ngal],
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t')

def build_corrected_randoms(catalog='LasDamasGeo', catalog_param_list={'n_mock':[[1, 'a']]}, corr='peak', corr_param={'sigma':5.3, 'fpeak':0.1},
    sanitycheck=False): 
    ''' 
    Construct randoms for fibercollision corrected data  
    '''
    # read series of corrected data in order to minimize single catalog dependent variance?
    # especially for las damas, which is divided into 4 parts 
    for i_catalog, catalog_param in enumerate(catalog_param_list['n_mock']): 
        corr_data = galaxy_data(DorR='data', catalog=catalog, catalog_param={'n_mock':catalog_param}, corr=corr, corr_param=corr_param) 
        print corr_data.file_name
        # target weights based on the corrected weights
        if i_catalog == 0: 
            if catalog.lower() == 'lasdamasgeo': 
                targ_weight = corr_data.weight 
                targ_z = corr_data.z
        else: 
            if catalog.lower() == 'lasdamasgeo': 
                targ_weight = np.concatenate((targ_weight, corr_data.weight))
                targ_z = np.concatenate((targ_z, corr_data.z))

    targ_weight_cum = targ_weight.cumsum()            # cumulative weights
    targ_weight_max = targ_weight_cum.max()           # max weight (a.k.a. total weight) 
    i_targ = np.arange(0,len(targ_weight_cum))
    print 'Corrected Galaxy Weight Sum = ', targ_weight_max

    # read true random
    true_random = galaxy_data(DorR='random', catalog=catalog, corr='true')
    Nran = len(true_random.ra)  # number of random galaxies
    print 'N_random = ', Nran

    corr_random_z = np.zeros(Nran) 
    # divide into chunks like make_catalog_z.py
    nchunk = 50 
    for i_chunk in range(nchunk): 
        # chunk indicies
        istart = Nran/(nchunk)*i_chunk
        iend = Nran/(nchunk)*(i_chunk+1)
        if i_chunk == nchunk-1: 
            iend = Nran 
        
        # randomly sample the sum of the target weights 
        wtarg = np.random.random(iend-istart)*targ_weight_max

        # get zindex of the target 
        zindx = np.floor(np.interp(wtarg, targ_weight_cum, i_targ)).astype(int)+1       # identical to make_catalog_z.py
        qqq = np.where(wtarg < targ_weight_cum[0])[0]
        zindx[qqq]  = 0 

        # assign redshift 
        corr_random_z[istart:iend] = targ_z[zindx]

    if len(corr_random_z) != Nran: 
        raise NameError("NOT ENOUGH CORRECTED RANDOMS!") 

    ## write Corrected Random File  
    # specify which correction scheme 
    if corr.lower() == 'upweight': 
        corr_str = '.upweight'
    #peak correction 
    elif (corr.lower() == 'peak') or (corr.lower() == 'peaknbar'): 
        corr_str = '.peak.sigma'+str(corr_param['sigma'])+'.fpeak'+str(corr_param['fpeak'])
    elif corr.lower() == 'peaktest':
        corr_str = '.peaktest.sigma'+str(corr_param['sigma'])+'.fpeak'+str(corr_param['fpeak'])

    # string to indicate all the corrected mocks files used to generate the randoms (e.g. 1a1b1c1d2a2b2c2d)
    if catalog.lower() == 'lasdamasgeo': 
        catalog_str = '.'
        if len(catalog_param_list['n_mock']) == 160:      # in the case that all mocks are used 
            catalog_str = '.allmocks'
        else: 
            for catalog_param in catalog_param_list['n_mock']: 
                catalog_str = ''.join([catalog_str, str(catalog_param[0]), catalog_param[1]])

    # corrected random file name 
    corr_random_filename = '/mount/riachuelo1/hahn/data/LasDamas/Geo/sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat'+\
            catalog_str+corr_str 
    # write! 
    np.savetxt(corr_random_filename, np.c_[true_random.ra, true_random.dec, corr_random_z], 
            fmt=['%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

def build_peakcorr_fibcoll(catalog='LasDamasGeo', catalog_param_list={'n_mock':[[1,'a']]}, corr='peaktest', corr_param={'sigma':5.3, 'fpeak':0.6}, 
        sanitycheck=False): 
    '''
    Apply peak correction to the fibercollided galaxy catalog 
    '''
    # cosmology is hardcoded based on fortran code cosmology 
    cosmo = {}
    cosmo['omega_M_0'] = 0.27 
    cosmo['omega_lambda_0'] = 0.73 
    cosmo['h'] = 0.7 
    cosmo = cosmos.distance.set_omega_k_0(cosmo) 

    if catalog.lower() == 'lasdamasgeo': 
        # LOS comoving distance for redshift limit
        comdis_lo = cosmos.distance.comoving_distance(0.16, **cosmo) 
        comdis_hi = cosmos.distance.comoving_distance(0.44, **cosmo) 

        # Loop through catalog parameter list  
        for n_catalog in catalog_param_list['n_mock']: 
            # read in fiber collided galaxies
            fibcoll_mock = galaxy_data(DorR='data', catalog='LasDamasGeo', catalog_param={'n_mock': n_catalog}, corr='Upweight') 
            appended_ra = [] 
            appended_dec = [] 
            appended_z = [] 
            appended_weight = [] 
            
            # if we want to check that the peak p(r) is generated properly
            if sanitycheck == True: 
                pr_test = [] 

            for i_mock in range(len(fibcoll_mock.weight)): 
                # for galaxies with wcp > 1
                while fibcoll_mock.weight[i_mock] > 1: 
                    fibcoll_mock.weight[i_mock] = fibcoll_mock.weight[i_mock]-1.0

                    # LOS comoving distance of the galaxy 
                    comdis_imock = cosmos.distance.comoving_distance(fibcoll_mock.z[i_mock], **cosmo) 
                    
                    rand_num = np.random.random(1) 
                    # if in the peak 
                    if rand_num < corr_param['fpeak']:          
                        # keep ra and dec
                        appended_ra.append(fibcoll_mock.ra[i_mock])
                        appended_dec.append(fibcoll_mock.dec[i_mock])

                        # appended galaxy has weight of 1.0 
                        appended_weight.append(1.0)
                    
                        # compute the displacement within peak ----------------------------------
                        rand1 = np.random.random(1) 
                        rand2 = np.random.random(1) 

                        rand2 = (-3.0+rand2*6.0)*corr_param['sigma']
                        peakpofr = np.exp(-1.0*np.abs(rand2)/corr_param['sigma'])
                        
                        while peakpofr <= rand1: 
                            rand1 = np.random.random(1) 
                            rand2 = np.random.random(1) 

                            rand2 = (-3.0+rand2*6.0)*corr_param['sigma']
                            peakpofr = np.exp(-1.0*np.abs(rand2)/corr_param['sigma'])
                        #----------------------------------- --------------------------------- 
                        # in case the displacement falls out of bound 
                        if (comdis_imock+rand2 > comdis_hi) or (comdis_imock+rand2 < comdis_lo): 
                            collided_z = comdis2z(comdis_imock-rand2)
                        else: 
                            collided_z = comdis2z(comdis_imock+rand2)

                        appended_z.append(collided_z[0]) 

                        if sanitycheck == True: 
                            pr_test.append(rand2) 
                    else:
                        if corr == 'peaktest': 
                            # do nothing 
                            pass
                        elif (corr == 'peak') or (corr == 'peaknbar'): 
                            # randomly displace 
                            appended_ra.append(fibcoll_mock.ra[i_mock])
                            appended_dec.append(fibcoll_mock.dec[i_mock])
                            appended_weight.append(1.0) 
                            
                            rand1 = np.random.random(1) 
                            appended_z.append(0.16+rand1[0]*0.28)
           
            print len(appended_ra), ' galaxies were peak corrected'
            fibcoll_mock.ra = np.concatenate([fibcoll_mock.ra, appended_ra])
            fibcoll_mock.dec = np.concatenate([fibcoll_mock.dec, appended_dec])
            fibcoll_mock.weight = np.concatenate([fibcoll_mock.weight, appended_weight])
            fibcoll_mock.z = np.concatenate([fibcoll_mock.z, appended_z])

            peakcorr = galaxy_data(DorR='data', catalog='LasDamasGeo', catalog_param={'n_mock': n_catalog}, 
                    corr=corr, corr_param=corr_param, readdata=False) 
            peakcorr_file = peakcorr.file_name 
            np.savetxt(peakcorr_file, np.c_[fibcoll_mock.ra, fibcoll_mock.dec, fibcoll_mock.z, fibcoll_mock.weight], 
                    fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5f'], delimiter='\t') 

def residual(arr1, arr2): 
    if len(arr1) != len(arr2): 
        raise TypeError("Input array lengths do not match.")
    else: 
        resid = np.array([arr1[i]/arr2[i] for i in range(len(arr1))])
    return resid

def comdis2z(comdis): 
    '''
    Given comoving distance and set cosmology, determine z 
    '''
    # set comoslogy 
    cosmo = {}
    cosmo['omega_M_0'] = 0.27 
    cosmo['omega_lambda_0'] = 0.73 
    cosmo['h'] = 0.7 
    cosmo = cosmos.distance.set_omega_k_0(cosmo) 

    z_arr = np.array([0.0+0.05*np.float(i) for i in range(21)]) 
    dm_arr = cosmos.distance.comoving_distance(z_arr, **cosmo) 

    # use numpy interpolate
    z = np.interp(comdis, dm_arr, z_arr) 
    return z 

def classify_triangles(k1, k2, k3, triangle='equilateral'): 
    '''
    Given k1, k2, k3, returns indices for (k1,k2,k3) that satify the specified triangle type  
    '''
    maxk = np.array([np.max([k1[i], k2[i], k3[i]]) for i in range(len(k1))])
    mink = np.array([np.min([k1[i], k2[i], k3[i]]) for i in range(len(k1))])
    if triangle == 'equilateral':       # only keep equilateral triangles
        triangle_index = (k1 == k2) & (k2 == k3) 
    elif triangle == 'acute':           # acute triangle
        triangle_index = (k1**2 + k2**2 > k3**2) & (k2**2 + k3**2 > k1**2) & (k3**2 + k1**2 > k2**2)
    elif triangle == 'obtuse':          # obtuse triangle
        triangle_index = (k1**2 + k2**2 < k3**2) | (k2**2 + k3**2 < k1**2) | (k3**2 + k1**2 < k2**2)
    elif triangle == 'extended':        # extended triangle   
        triangle_index = maxk/mink > 3.0
    return triangle_index

def append_corr_nbar(DorR='data', catalog='LasDamasGeo', catalog_param={'n_mock':[1,'a']},
        corr='peaknbar', corr_param={'sigma':5.3, 'fpeak':0.5}, sanitycheck=False):
    '''
    append corrected interpolated nbar(z) to corrected data or random file
    '''
    # read in data/random galaxies
    gal_data = galaxy_data(DorR=DorR, catalog=catalog, catalog_param=catalog_param, corr=corr, corr_param=corr_param)
    if DorR == 'random':
        print 'Nran = ', len(gal_data.ra)
    # read in corrected nbar file
    corr_nbar = nbar(catalog=catalog, corr=corr, corr_param=corr_param)

    # interpolate within redshift limits
    if catalog.lower() == 'lasdamasgeo':
        zlim = (corr_nbar.zmid > 0.16) & (corr_nbar.zmid < 0.44)        # for las damas geo
    else:
        raise NameError("not yet coded!")
    # numpy interpolate
    nbar_arr = np.interp(gal_data.z, corr_nbar.zmid[zlim], corr_nbar.nbar[zlim])

    if DorR == 'random':
        print 'Nran = ', len(nbar_arr)

    # if santiy check is true then plot the interpolated nbar values
    if sanitycheck == True:
        fig = plt.figure(1)
        sub = fig.add_subplot(111)
        sub.scatter(gal_data.z, nbar_arr, color='red', s=5, label=DorR)
        sub.plot(corr_nbar.zmid, corr_nbar.nbar, color='black', lw=2, label="Corrected nbar(z)")
        sub.set_xlim([0.16, 0.44])
        sub.set_ylim([9.3e-05, 10.0e-05])
        sub.legend(loc='upper right', scatterpoints=1)
        fig.savefig(''.join(['/home/users/hahn/research/figures/boss/',
        'corrected-', DorR, '-nbarz-sanitycheck.png']), bbox_inches="tight")
        fig.clear()

    # write corr nbar appended data/random
    # corr nbar appended file name
    gal_corr_nbar_file = ''.join([gal_data.file_name, '.corrnbar'])
    if DorR == 'data':
        np.savetxt(gal_corr_nbar_file,
                np.c_[gal_data.ra, gal_data.dec, gal_data.z, nbar_arr, gal_data.weight],
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5e', '%10.5f'], delimiter='\t')
    elif DorR == 'random':
        np.savetxt(gal_corr_nbar_file,
                np.c_[gal_data.ra, gal_data.dec, gal_data.z, nbar_arr],
                fmt=['%10.5f', '%10.5f', '%10.5f', '%10.5e'], delimiter='\t')

def get_fibcoll_dir(file_type, catalog='lasdamasgeo'): 
    '''
    get data/FFT/power directories given catalog
    '''
    if (file_type.lower() != 'data') and (file_type.lower() != 'fft') and (file_type.lower() != 'power'): 
        raise NameError('either data, fft, or power') 
    else: 
        # for lasdamasgeo 
        if catalog.lower() == 'lasdamasgeo': 
            # lasdamasgeo data
            if file_type.lower() == 'data': 
                file_dir = '/mount/riachuelo1/hahn/data/LasDamas/Geo/'
            elif file_type.lower() == 'fft': 
                file_dir = '/mount/riachuelo1/hahn/FFT/LasDamas/Geo/'
            else:
                file_dir = '/mount/riachuelo1/hahn/power/LasDamas/Geo/'
    return file_dir 

# Plotting -----------------------------------------------------------------
def plot_powerspec_fibcolcorr_comparison(catalog='PTHalo', filespec={'version':'v11p0', 'n_mock':range(1,11), 'Nrandom':100}, 
        correction_method=['true', 'delta', 'peak'], corrspec={'sigma':[0, 0, 6.3], 'fpeak':[1.0, 1.0, 1.0]}, resid='False'):
    '''
    Comparison of P(k)avg for different fibercollision correciton methods
    Make sure 'true' is first in the correction method list
    '''
    prettyplot()
    
    fig = plt.figure(1, figsize=(7, 8))
    sub = fig.add_subplot(111)

    # Tableau20 colors 
    corr_color = [(89, 89, 89), (31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),  
            (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),  
            (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),  
            (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),  
            (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  
  
    # Scale the RGB values to the [0, 1] range
    for i in range(len(corr_color)):  
        r, g, b = corr_color[i]  
        corr_color[i] = (r / 255., g / 255., b / 255.)  
             
    sigma_flag = ''
    for i_corr, correction in enumerate(correction_method):     # loop through correction methods
        for i_mock in filespec['n_mock']:                       # compute average[P(k)] for each correction method
            i_filespec = filespec.copy()
            i_filespec['n_mock'] = i_mock 
            i_corrspec = corrspec.copy() 

            for key in corrspec.keys(): 
                i_corrspec[key] = corrspec[key][i_corr]
            power_i = spec(spectrum='power', catalog=catalog, filespec=i_filespec, fibcolcorr=correction, corrspec=i_corrspec) 
            power_i.readfile()
            if i_mock == filespec['n_mock'][0]: 
                avg_k = power_i.k
                sum_Pk = power_i.Pk
            else: 
                sum_Pk = sum_Pk + power_i.Pk
        avg_Pk = [ sum_Pk[i]/float(len(filespec['n_mock'])) for i in range(len(sum_Pk))]    # average P(k)
        
        # P(k) Comparison
        if resid == 'False':  
            if correction == 'true':
                sub.plot(avg_k, avg_Pk, color=corr_color[i_corr], label=r"$\mathtt{\overline{P(k)_{"+correction+"}}}$")
            elif correction == 'peak':
                sub.scatter(avg_k, avg_Pk, color=corr_color[i_corr], label=r"$\mathtt{\overline{P(k)_{"+correction+','+\
                        ''.join([str(corrspec['sigma'][i_corr]), str(corrspec['fpeak'][i_corr])])+"}}}$")
            else: 
                sub.scatter(avg_k, avg_Pk, color=corr_color[i_corr], label=r"$\mathtt{\overline{P(k)_{"+correction+"}}}$")

        # P(k) residual comparison
        else:               
            if correction == 'true':
                avg_Pk_true = avg_Pk
            elif correction == 'peak': 
                resid_label = ''.join([r"$\mathtt{", "\overline{P(k)_{", correction, ",", 
                    ''.join([str(corrspec['sigma'][i_corr]), str(corrspec['fpeak'][i_corr])]), "}}/\overline{P(k)_{True}}}$"])
                sub.scatter(avg_k, residual(avg_Pk, avg_Pk_true), color=corr_color[i_corr], label=resid_label)
            else:
                resid_label = ''.join([r"$\mathtt{", "\overline{P(k)_{",correction, "}}/\overline{P(k)_{True}}}$"])
                sub.scatter(avg_k, residual(avg_Pk, avg_Pk_true), color=corr_color[i_corr], label=resid_label)
        
        # sigma label 
        if correction == 'peak': 
            sigma_flag = ''.join([sigma_flag, '_sigma', str(corrspec['sigma'][i_corr]), '_fpeak', str(corrspec['fpeak'][i_corr])]) 

    if resid == 'False':
        if catalog.lower() == 'lasdamasgeo':
            ylimit = [10**2,10**5.5]
            ytext = 10**5.2
        elif catalog.lower() == 'tilingmock':
            ylimit = [10**2,10**5.5]
            ytext = 10**5.2
        elif catalog.lower() == 'pthalo':
            ylimit = [10**2,10**5.5]
            ytext = 10**5.2
        ylabel = 'P(k)'
        sub.set_yscale('log')
        fig_name = ''.join(['powerspec_', catalog.lower(), 
            '_fibcoll_', '_'.join(correction_method), sigma_flag, '_comparison.png'])
    else:
        ylimit = [0.9,1.2]
        ytext = 1.15
        ylabel = r"$\mathtt{\overline{P(k)}/\overline{P(k)_{\rm{True}}}}$"
        fig_name = ''.join(['powerspec_', catalog.lower(), 
            '_fibcoll_', '_'.join(correction_method), sigma_flag, '_residual_comparison.png'])     

    sub.text(1.5*10**-3, np.mean(ylimit), 
            ''.join([str(len(filespec['n_mock'])), ' ', catalog.upper()]))              # number of mocks + Catalog name 
    sub.set_xscale('log')
    sub.set_xlim([10**-3,10**0])
    sub.set_ylim(ylimit)
    sub.set_xlabel('k', fontsize=20)
    sub.set_ylabel(ylabel, fontsize=20)

    if resid == 'False': 
        sub.legend(loc='lower left', scatterpoints=1, prop={'size':14})
    elif resid == 'True': 
        sub.legend(loc='upper left', scatterpoints=1, prop={'size':14})

    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', fig_name]), bbox_inches="tight")
    fig.clear()

def plot_bispec_fibcolcorr_comparison(BorQ='B', x_axis='triangles', triangle='all', 
        catalog='PTHalo', filespec={'version':'v11p0', 'n_mock':range(1,11), 'Nrandom':100}, 
        correction_method=['true', 'delta', 'peak'], corrspec={'sigma':[0, 0, 6.3], 'fpeak':[1.0, 1.0, 1.0]}, resid='False'):
    '''
    Comparison of B(k)avg for different fibercollision correciton methods
    Make sure 'true' is first in the correction method list
    '''
    prettyplot()
    
    fig = plt.figure(1, figsize=(7, 8))
    sub = fig.add_subplot(111)

    corr_color = ['grey', 'blue', 'red']
    for i_corr, correction in enumerate(correction_method):     # loop through correction methods
        for i_mock in filespec['n_mock']:                       # compute average[B(k)] for each correction method
            i_filespec = filespec.copy()
            i_filespec['n_mock'] = i_mock 
            i_corrspec = corrspec.copy() 
            for key in corrspec.keys(): 
                i_corrspec[key] = corrspec[key][i_corr]

            bispec_i = spec(spectrum='bispec', catalog=catalog, filespec=i_filespec, fibcolcorr=correction, corrspec=i_corrspec) 
            bispec_i.readfile()
            if i_mock == filespec['n_mock'][0]: 
                if x_axis == 'triangles':                       # specify x-axis (triangle index, avgk, or max k) 
                    avg_k = bispec_i.i_triangle                 # x-axis is triangles
                elif x_axis == 'avg_k': 
                    avg_k = bispec_i.avgk                       # x-axis is avg(k1,k2,k3)
                elif x_axis == 'max_k': 
                    avg_k = bispec_i.kmax                       # x-axis is max(k1,k2,k3) 
                avg_k1 = bispec_i.k1
                avg_k2 = bispec_i.k2
                avg_k3 = bispec_i.k3

                if BorQ == 'B': 
                    sum_Bk = bispec_i.Bk
                elif BorQ == 'Q': 
                    sum_Bk = bispec_i.Q
            else: 
                if BorQ == 'B': 
                    sum_Bk = sum_Bk + bispec_i.Bk
                elif BorQ == 'Q': 
                    sum_Bk = sum_Bk + bispec_i.Q
        avg_Bk = [ sum_Bk[i]/float(len(filespec['n_mock'])) for i in range(len(sum_Bk))]    # average B(k)

        if triangle == 'all': 
            tri_index = np.array([True for i in range(len(avg_Bk))], dtype=bool)
        else: 
            tri_index = classify_triangles(avg_k1, avg_k2, avg_k3, triangle=triangle)
            print tri_index

        if resid == 'False':                                    # B(k) Comparison
            if correction == 'true':
                sub.scatter(avg_k[tri_index], avg_Bk[tri_index], 
                        color=corr_color[i_corr], label=r"$\mathtt{\overline{"+BorQ+"(k)_{"+correction+"}}}$")
            else:
                sub.scatter(avg_k[tri_index], avg_Bk[tri_index], 
                        color=corr_color[i_corr], label=r"$\mathtt{\overline{"+BorQ+"(k)_{"+correction+"}}}$")
        else:                                                   # B(k) residual comparison
            if correction == 'true':
                avg_Bk_true = avg_Bk
            else:
                sub.scatter(np.array(avg_k)[tri_index], residual(np.array(avg_Bk)[tri_index], np.array(avg_Bk_true)[tri_index]), 
                        color=corr_color[i_corr], 
                        label=''.join([r"$\mathtt{", "\overline{", BorQ, "(k)_{",correction, "}}/\overline{", BorQ, "(k)_{True}}}$"]))
        if correction == 'peak': 
            sigma_flag = ''.join(['_sigma', str(corrspec['sigma'][i_corr]), '_fpeak', str(corrspec['fpeak'][i_corr])]) 
        else: 
            sigma_flag = ''

    if resid=='False':
        if catalog.lower() == 'tilingmock':
            ylimit = [10**2,10**5.5]
            ytext = 10**5.2
        elif catalog.lower() == 'pthalo':
            ylimit = [10**2,10**5.5]
            ytext = 10**5.2
        ylabel = '$\mathtt{'+BorQ+'(k)}$'
        sub.set_yscale('log')
        fig_name = ''.join(['bispec_', BorQ, 'k_', x_axis, '_', catalog.lower(), 
            '_fibcoll_', '_'.join(correction_method), sigma_flag, '_', triangle, 'triangles_comparison.png'])
    else:
        ylimit = [0.9,1.2]
        ytext = 1.15
        ylabel = r"$\mathtt{\overline{"+BorQ+"(k)}/\overline{"+BorQ+"(k)_{\rm{True}}}}$"
        fig_name = ''.join(['bispec_', BorQ, 'k_', x_axis, '_', catalog.lower(), 
            '_fibcoll_', '_'.join(correction_method), sigma_flag, '_', triangle, 'triangles_residual_comparison.png'])     
        sub.set_ylim(ylimit)
    sub.text(1.5*10**-3,ytext, 
            ''.join([str(len(filespec['n_mock'])), ' ', catalog.upper()]))              # number of mocks + Catalog name 
    if x_axis == 'triangles': 
        sub.set_xlabel('Triangles', fontsize=20)
        sub.set_xlim([0, 7500])
    elif x_axis == 'avg_k': 
        sub.set_xlabel('avg(k1, k2, k3)', fontsize=20)
        sub.set_xscale('log')
        sub.set_xlim([10**-3,10**0])
    elif x_axis == 'max_k': 
        sub.set_xlabel('max(k1, k2, k3)', fontsize=20)
        sub.set_xscale('log')
        sub.set_xlim([10**-3,10**0])
    sub.set_ylabel(ylabel, fontsize=20)
    sub.legend(loc='lower left', prop={'size':14}, scatterpoints=1)
    sub.grid(True)
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/fiber_collision/', fig_name]))
    fig.clear()

def plot_corrected_random_nbar_comp(catalog='lasdamasgeo', catalog_param_list=[[1, 'a']], corr='peak', corr_param=[5.3, 0.1]): 
    '''
    Sanity check to check whether or not the nbar(z) for the corrected randoms are adjusted in the same way as the corrected
    mocks
    '''
    # import arbitrary nbar file for z values
    zcen, zlow, zhigh = np.loadtxt("/mount/riachuelo1/hahn/data/nbar-junk.dat" , unpack=True, usecols=[0,1,2])
    
    # configure plot
    prettyplot()
    fig = plt.figure(1, figsize=(7,8))
    sub = fig.add_subplot(111)
    
    # read true random
    true_random = galaxy_data(DorR='random', catalog=catalog, corr='true')
    Nran = np.float(len(true_random.ra))  # number of random galaxies
    print true_random.filename
    print 'N_random = ', Nran

    # read corrected random 
    corr_random = galaxy_data(DorR='random', catalog=catalog, corr='peak')
    Nran_corr = np.float(len(corr_random.ra))
    print corr_random.filename
    print 'N_random,corr = ', Nran_corr
   
    # calculate the random file's effective nbar
    n_rand_true = np.zeros(len(zcen))
    n_rand_corr = np.zeros(len(zcen))
    for i_z, zmid in enumerate(zcen): 
        true_rand_z_lim = (true_random.z >= zlow[i_z]) & (true_random.z < zhigh[i_z])
        n_rand_true[i_z] = np.float(len(true_random.z[true_rand_z_lim])) 
        
        corr_rand_z_lim = (corr_random.z >= zlow[i_z]) & (corr_random.z < zhigh[i_z])
        n_rand_corr[i_z] = np.float(len(corr_random.z[corr_rand_z_lim]))
    
    true_nozero = n_rand_true > 0.0
    
    # read series of corrected data
    for i_catalog, catalog_param in enumerate(catalog_param_list): 
        corr_data = galaxy_data(DorR='data', catalog=catalog, catalog_param=catalog_param, corr=corr, corr_param=corr_param) 
        corr_data_wsum = np.sum(corr_data.weight)
        print corr_data.filename
        print corr_data_wsum
        
        n_gal_corr = np.zeros(len(zcen))
        n_gal_corr_scaled = np.zeros(len(zcen))
        for i_z, zmid in enumerate(zcen): 
            corr_z_lim = (corr_data.z >= zlow[i_z]) & (corr_data.z < zhigh[i_z])
            n_gal_corr[i_z] = np.sum(corr_data.weight[corr_z_lim])
            n_gal_corr_scaled[i_z] = (n_gal_corr[i_z]/corr_data_wsum)*np.float(Nran)

        sub.plot(zcen[true_nozero], n_gal_corr_scaled[true_nozero]/n_rand_true[true_nozero], 
                color='grey', lw=3, label=r'$corr_{data}/random$') 
        sub.plot(zcen[true_nozero], n_gal_corr_scaled[true_nozero]/n_rand_corr[true_nozero], 
                color='blue', lw=3, label=r'$corr_{data}/corr_{random}$') 

    sub.plot(zcen[true_nozero], n_rand_corr[true_nozero]/n_rand_true[true_nozero], color='red', lw=3, label=r'$corr_{random}/random$') 

    sub.set_ylabel(r'$\bar{n}(z)$ ratio')
    sub.set_xlim([0.16, 0.44])
    sub.set_ylim([0.9, 1.3])
    sub.set_xlabel('z')
    sub.legend(loc='best')
    fig.savefig(''.join(['/home/users/hahn/research/figures/boss/', 
        'lasdamasgeo_corrected_random_nbar_comparison.png']), bbox_inches="tight")
    fig.clear() 

# Runs -----------------------------------------------------------------
def lasdamasgeo_fibcoll_pk(n_mocks, corr='peak'): 
    '''
    Compute fibercollision corrected P(k)
    '''
    # some constants throughout the code
    P0 = 20000
    sscale=3600.0
    Rbox=1800.0
    box="3600"
    grid="360"
    sigma = 5.3
    peak_fracs = [0.1] #[0.1*np.float(i) for i in range(11)]
    nbar_file = "/mount/riachuelo1/hahn/data/nbar-junk.dat"             # just some junk nbar file 

    # file part strings 
    mock_prefix = 'sdssmock_gamma_lrgFull_zm_oriana'
    mock_suffix = '_no.rdcz.dat'
    fibcoll_suffix = '_no.rdcz.fibcoll.dat'

    # compile all Fortran codes 
    ldg_code_dir = '/home/users/hahn/powercode/FiberCollision/LasDamas/Geo/'
    peakcorr_fortcode = 'peakcorr_ldg_fibcollmock.f'
    fft_w_360grid_fortcode = 'FFT_ldg_fkp_w_360grid.f'
    power_360grid_fortcode = 'power_ldg_fkp_360grid_180bin.f'

    # for each fortran code
    print "Compiling Fortran Codes"
    for code in [peakcorr_fortcode, fft_w_360grid_fortcode, power_360grid_fortcode]: 
        # executable file 
        fort_exe = 'exe/'+code.rsplit('.')[0]+'.exe'
        # compile command
        compile_cmd = ' '.join(['ifort -O3 -o', ldg_code_dir+fort_exe, ldg_code_dir+code, 
            '-L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw'])
        print compile_cmd
        # call compile command 
        subprocess.call(compile_cmd.split())
   
    # For upweight or peak correction method 
    if (corr.lower() == 'upweight') or (corr.lower() == 'peak') or (corr.lower() == 'peaknbar') or (corr.lower() == 'peaktest'): 
        # generate file lists to call them easily later  
        file_indices = []
        fibercollided_files = []
        for i_mock in range(1, n_mocks+1): 
            for letter in ['a', 'b', 'c', 'd']: 
                # append file indicies 
                file_indices.append([i_mock, letter])
    
                # append fibercollided file names 
                fibercollided = galaxy_data(DorR='data', catalog='LasDamasGeo', catalog_param={'n_mock': [i_mock, letter]}, 
                        corr='upweight', readdata=False) 
                fibercollided_files.append(fibercollided.file_name)

        for i_file, fibercollided_file in enumerate(fibercollided_files):  
            # fibercollided file does NOT exist
            if os.path.isfile(fibercollided_file) == False:  
                fibcollided_cmd = ' '.join(['idl', '-e', ldg_code_dir+'ldg_fibcollmock_wcp_assign,', 
                    str(file_indices[i_file][0]), ",'"+file_indices[i_file][1]+"'"])
                print "Fibercolliding the Mocks!"
                subprocess.call(fibcollided_cmd.split())
            else: 
                print fibercollided_file, ' already exists'
        
        # For upweight fiber-collision correction method 
        if corr.lower() == 'upweight': 
            # append corrected nbar to peak corrected mocks
            for i_file, file_index in enumerate(file_indices): 
                if os.path.isfile(fibercollided_files[i_file]+'.corrnbar') == False:        # if there is no corrected mock file
                    print "appending corrected nbar to mock!"
                    append_corr_nbar(DorR='data', catalog='LasDamasGeo', catalog_param={'n_mock':[file_index[0],file_index[1]]},
                            corr='upweight', sanitycheck=False)
                    fibercollided_files[i_file] = fibercollided_files[i_file]+'.corrnbar'
                else: 
                    fibercollided_files[i_file] = fibercollided_files[i_file]+'.corrnbar'
                    print fibercollided_files[i_file], ' already exists'

            # Corrected RANDOMS --------------------------------------------------------------------------------------------
            corr_rand_file = ''.join([get_fibcoll_dir('data', catalog='lasdamasgeo'), 
                'sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat.allmocks.upweight']) 
            if os.path.isfile(corr_rand_file) == False:                         # if there is no corrected random file
                print "Building Correct Randoms! (Will take a while!)" 
                build_ldg_upweight_corrected_randoms() 
            else: 
                print corr_rand_file, " already exists!"

            # append corrected nbar to corrected random 
            if os.path.isfile(corr_rand_file+'.corrnbar') == False: 
                print "appending corrected nbar to random!"
                append_corr_nbar(DorR='random', 
                        catalog='LasDamasGeo', catalog_param={'n_mock':[1,'a']},
                        corr='upweight', sanitycheck=False)
                corr_rand_file = corr_rand_file+'.corrnbar'
            else: 
                corr_rand_file = corr_rand_file+'.corrnbar'
                print corr_rand_file, " already exists!"

            # FFT ---------------------------------------------------------------------------------------------------------
            FFT_exe = ldg_code_dir+'exe/'+fft_w_360grid_fortcode.rsplit('.')[0]+'.exe'

            # FFT for random 
            fft_rand_file = ''.join([get_fibcoll_dir('fft', catalog='lasdamasgeo'), 
                'FFT_', corr_rand_file.rsplit('/')[-1], '.grid', grid, '.P0', str(P0), '.box', box])
            # FFT random command 
            FFT_rand_cmd = ' '.join([FFT_exe, str(Rbox), "1", str(P0), corr_rand_file, fft_rand_file])
            # call FFT randomc ommand 
            if os.path.isfile(fft_rand_file) == False: 
                print "Building ", fft_rand_file 
                subprocess.call(FFT_rand_cmd.split())
            else: 
                print fft_rand_file, " already exists" 
            
            # FFT for mocks
            fft_files = [] 
            for i_file, file_index in enumerate(file_indices): 
                fft_file = ''.join([get_fibcoll_dir('fft', catalog='lasdamasgeo'), 
                    'FFT_', fibercollided_files[i_file].rsplit('/')[-1], '.upweight.grid', grid, '.P0', str(P0), '.box', box])
                fft_files.append(fft_file) 
                
                # FFT command
                print 'Building ', fft_file 
                FFT_cmd = ' '.join([FFT_exe, str(Rbox), "0", str(P0), fibercollided_files[i_file], fft_file]) 
                subprocess.call(FFT_cmd.split()) 

            # P(k) ---------------------------------------------------------------------------------------------------------
            power_exe = ldg_code_dir+'exe/'+power_360grid_fortcode.rsplit('.')[0]+'.exe'
            for i_file, file_index in enumerate(file_indices): 
                power_file = ''.join([get_fibcoll_dir('power', catalog='lasdamasgeo'), 
                    'power_', fibercollided_files[i_file].rsplit('/')[-1], '.upweight.grid', grid, '.P0', str(P0), '.box', box])
                
                power_cmd = ' '.join([power_exe, fft_files[i_file], fft_rand_file, power_file, str(sscale)]) 
                print 'Building ', power_file 
                subprocess.call(power_cmd.split())
        
        elif (corr.lower() == 'peak') or (corr.lower() == 'peaknbar'): 
            # loop for peak fractions 
            for peak_frac in peak_fracs: 
                # correction string
                corr_str = '.peak.sigma'+str(sigma)+'.fpeak'+str(peak_frac)

                # First MOCKS --------------------------------------------------------------------------------------------
                # correct the fibercollided mocks using peak correction method
                peakcorrected_files = [] 
                for i_file, file_index in enumerate(file_indices): 
                    # store corrected file names 
                    peakcorrected_file = ''.join([get_fibcoll_dir('data', catalog='lasdamasgeo'),
                        mock_prefix, str(100+file_index[0])[1:3], file_index[1], fibcoll_suffix, corr_str])
                    peakcorrected_files.append(peakcorrected_file)

                    if os.path.isfile(peakcorrected_file) == False:                     # if there is no corrected mock file
                        corr_cmd = ' '.join([ldg_code_dir+'exe/'+peakcorr_fortcode.rsplit('.')[0]+'.exe', 
                            str(sigma), str(peak_frac), nbar_file, fibercollided_files[i_file], peakcorrected_file])
                        print "Peak Correcting the Fiber collided mocks with sigma = ", sigma, " and fpeak = ", peak_frac
                        subprocess.call(corr_cmd.split())
                    else: 
                        print peakcorrected_file, ' already exists'

                # append corrected nbar to peak corrected mocks
                for i_file, file_index in enumerate(file_indices): 
                    if os.path.isfile(peakcorrected_files[i_file]+'.corrnbar') == False:        # if there is no corrected mock file
                        print "appending corrected nbar to mock!"
                        append_corr_nbar(DorR='data', catalog='LasDamasGeo', catalog_param={'n_mock':[file_index[0],file_index[1]]},
                                corr='peaknbar', corr_param={'sigma':sigma, 'fpeak':peak_frac}, sanitycheck=False)
                        peakcorrected_files[i_file] = peakcorrected_files[i_file]+'.corrnbar'
                    else: 
                        peakcorrected_files[i_file] = peakcorrected_files[i_file]+'.corrnbar'
                        print peakcorrected_files[i_file], ' already exists'

                # Then RANDOMS --------------------------------------------------------------------------------------------
                corr_rand_file = ''.join([get_fibcoll_dir('data', catalog='lasdamasgeo'), 
                    'sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat.allmocks.peak.sigma', str(sigma), '.fpeak', str(peak_frac)]) 
                if os.path.isfile(corr_rand_file) == False:                         # if there is no corrected random file
                    print "Building Correct Randoms! (Will take a while!)" 
                    build_ldg_peak_corrected_randoms(sigma, peak_frac, corr='peak') 
                else: 
                    print corr_rand_file, " already exists!"

                # append corrected nbar to corrected random 
                if os.path.isfile(corr_rand_file+'.corrnbar') == False: 
                    print "appending corrected nbar to random!"
                    append_corr_nbar(DorR='random', 
                            catalog='LasDamasGeo', catalog_param={'n_mock':[1,'a']},
                            corr='peaknbar', corr_param={'sigma':sigma, 'fpeak':peak_frac}, 
                            sanitycheck=False)
                    corr_rand_file = corr_rand_file+'.corrnbar'
                else: 
                    corr_rand_file = corr_rand_file+'.corrnbar'
                    print corr_rand_file, " already exists!"

                # FFT ---------------------------------------------------------------------------------------------------------
                FFT_exe = ldg_code_dir+'exe/'+fft_w_360grid_fortcode.rsplit('.')[0]+'.exe'

                # FFT for random 
                fft_rand_file = ''.join([get_fibcoll_dir('fft', catalog='lasdamasgeo'), 
                    'FFT_', corr_rand_file.rsplit('/')[-1], '.grid', grid, '.P0', str(P0), '.box', box])
                FFT_rand_cmd = ' '.join([FFT_exe, str(Rbox), "1", str(P0), corr_rand_file, fft_rand_file])
                if os.path.isfile(fft_rand_file) == False: 
                    print "Building ", fft_rand_file 
                    subprocess.call(FFT_rand_cmd.split())
                else: 
                    print fft_rand_file, " already exists" 
                
                # FFT for mocks
                fft_files = [] 
                for i_file, file_index in enumerate(file_indices): 
                    fft_file = ''.join([get_fibcoll_dir('fft', catalog='lasdamasgeo'), 
                        'FFT_', peakcorrected_files[i_file].rsplit('/')[-1], '.grid', grid, '.P0', str(P0), '.box', box])
                    fft_files.append(fft_file) 
                    
                    # FFT command
                    print 'Building ', fft_file 
                    FFT_cmd = ' '.join([FFT_exe, str(Rbox), "0", str(P0), peakcorrected_files[i_file], fft_file]) 
                    subprocess.call(FFT_cmd.split()) 

                # P(k) ---------------------------------------------------------------------------------------------------------
                power_exe = ldg_code_dir+'exe/'+power_360grid_fortcode.rsplit('.')[0]+'.exe'
                for i_file, file_index in enumerate(file_indices): 
                    power_file = ''.join([get_fibcoll_dir('power', catalog='lasdamasgeo'), 
                        'power_', peakcorrected_files[i_file].rsplit('/')[-1], '.grid', grid, '.P0', str(P0), '.box', box])
                    
                    power_cmd = ' '.join([power_exe, fft_files[i_file], fft_rand_file, power_file, str(sscale)]) 
                    print 'Building ', power_file 
                    subprocess.call(power_cmd.split())

        elif corr.lower() == 'peaktest': 
            print 'PEAK TEST CORRECTION'
            # loop for peak fractions 
            for peak_frac in peak_fracs: 
                # correction string
                corr_str = '.peaktest.sigma'+str(sigma)+'.fpeak'+str(peak_frac)

                # First MOCKS --------------------------------------------------------------------------------------------
                # correct the fibercollided mocks using peak correction method
                peakcorrected_files = [] 
                for i_file, file_index in enumerate(file_indices): 
                    # store corrected file names 
                    peakcorrected = galaxy_data(DorR='data', catalog='LasDamasGeo', 
                            catalog_param={'n_mock': [file_index[0], file_index[1]]}, 
                            corr='peaktest', corr_param={'sigma': sigma, 'fpeak':peak_frac}, readdata=False) 
                    peakcorrected_file = peakcorrected.file_name
                    peakcorrected_files.append(peakcorrected_file) 

                    if os.path.isfile(peakcorrected_file) == False:                     # if there is no corrected mock file
                        print "Peak Test Correcting the Fiber collided mocks with sigma = ", sigma, " and fpeak = ", peak_frac
                        build_peakcorr_fibcoll(catalog='LasDamasGeo', catalog_param_list={'n_mock':[[file_index[0], file_index[1]]]}, 
                            corr_param={'sigma':sigma, 'fpeak':peak_frac}) 
                    else: 
                        print peakcorrected_file, ' already exists'

                # append corrected nbar to peak corrected mocks
                for i_file, file_index in enumerate(file_indices): 
                    if os.path.isfile(peakcorrected_files[i_file]+'.corrnbar') == False:        # if there is no corrected mock file
                        print "appending corrected nbar to mock!"
                        append_corr_nbar(DorR='data', catalog='LasDamasGeo', catalog_param={'n_mock':[file_index[0],file_index[1]]},
                                corr='peaktest', corr_param={'sigma':sigma, 'fpeak':peak_frac}, sanitycheck=False)
                        peakcorrected_files[i_file] = peakcorrected_files[i_file]+'.corrnbar'
                    else: 
                        peakcorrected_files[i_file] = peakcorrected_files[i_file]+'.corrnbar'
                        print peakcorrected_files[i_file], ' already exists'

                # Then RANDOMS --------------------------------------------------------------------------------------------
                corr_rand = galaxy_data(DorR='random', catalog='LasDamasGeo', 
                        corr='peaktest', corr_param={'sigma': sigma, 'fpeak':peak_frac}, readdata=False) 
                corr_rand_file = corr_rand.file_name
                if os.path.isfile(corr_rand_file) == False:                         # if there is no corrected random file
                    print "Building Correct Randoms! (Will take a while!)" 
                    build_ldg_peak_corrected_randoms(sigma, peak_frac, corr='peaktest') 
                else: 
                    print corr_rand_file, " already exists!"

                # append corrected nbar to corrected random 
                if os.path.isfile(corr_rand_file+'.corrnbar') == False: 
                    print "appending corrected nbar to random!"
                    append_corr_nbar(DorR='random', 
                            catalog='LasDamasGeo', catalog_param={'n_mock':[1,'a']},
                            corr='peaktest', corr_param={'sigma':sigma, 'fpeak':peak_frac}, 
                            sanitycheck=False)
                    corr_rand_file = corr_rand_file+'.corrnbar'
                else: 
                    corr_rand_file = corr_rand_file+'.corrnbar'
                    print corr_rand_file, " already exists!"

                # FFT ---------------------------------------------------------------------------------------------------------
                FFT_exe = ldg_code_dir+'exe/'+fft_w_360grid_fortcode.rsplit('.')[0]+'.exe'

                # FFT for random 
                fft_rand_file = ''.join([get_fibcoll_dir('fft', catalog='lasdamasgeo'), 
                    'FFT_', corr_rand_file.rsplit('/')[-1], '.grid', grid, '.P0', str(P0), '.box', box])
                FFT_rand_cmd = ' '.join([FFT_exe, str(Rbox), "1", str(P0), corr_rand_file, fft_rand_file])
                if os.path.isfile(fft_rand_file) == False: 
                    print "Building ", fft_rand_file 
                    subprocess.call(FFT_rand_cmd.split())
                else: 
                    print fft_rand_file, " already exists" 
                
                # FFT for mocks
                fft_files = [] 
                for i_file, file_index in enumerate(file_indices): 
                    fft_file = ''.join([get_fibcoll_dir('fft', catalog='lasdamasgeo'), 
                        'FFT_', peakcorrected_files[i_file].rsplit('/')[-1], '.grid', grid, '.P0', str(P0), '.box', box])
                    fft_files.append(fft_file) 
                    
                    # FFT command
                    print 'Building ', fft_file 
                    FFT_cmd = ' '.join([FFT_exe, str(Rbox), "0", str(P0), peakcorrected_files[i_file], fft_file]) 
                    subprocess.call(FFT_cmd.split()) 

                # P(k) ---------------------------------------------------------------------------------------------------------
                power_exe = ldg_code_dir+'exe/'+power_360grid_fortcode.rsplit('.')[0]+'.exe'
                for i_file, file_index in enumerate(file_indices): 
                    power_file = ''.join([get_fibcoll_dir('power', catalog='lasdamasgeo'), 
                        'power_', peakcorrected_files[i_file].rsplit('/')[-1], '.grid', grid, '.P0', str(P0), '.box', box])
                    
                    power_cmd = ' '.join([power_exe, fft_files[i_file], fft_rand_file, power_file, str(sscale)]) 
                    print 'Building ', power_file 
                    subprocess.call(power_cmd.split())

def lasdamasgeo_fibcoll_pk_rand_disp_test(n_mocks): 
    '''
    Random displacement test for lasdamasgeo P(k)
    '''
    # some constants throughout the code
    P0 = 20000
    sscale=3600.0
    Rbox=1800.0
    box="3600"
    grid="360"
    sigma = 5.3
    nbar_file = "/mount/riachuelo1/hahn/data/nbar-junk.dat"             # just some junk nbar file 

    # file part strings 
    mock_prefix = 'sdssmock_gamma_lrgFull_zm_oriana'
    mock_suffix = '_no.rdcz.dat'
    test_suffix = '_no.rdcz.randdisptest.dat'

    rand_file = '/mount/chichipio2/rs123/MOCKS/randoms/sdssmock_gamma_lrgFull.rand_200x_no.rdcz.dat'

    # Fortran code compilation -------------------------------------------------------------------- 
    # compile all Fortran codes 
    ldg_code_dir = '/home/users/hahn/powercode/FiberCollision/LasDamas/Geo/'
    peakcorr_fortcode = 'randdisp_ldg_test.f'
    fft_w_360grid_fortcode = 'FFT_ldg_fkp_rand_disp_test_360grid.f'
    power_360grid_fortcode = 'power_ldg_fkp_360grid_180bin.f'

    # for each fortran code
    print "Compiling Fortran Codes"
    for code in [peakcorr_fortcode, fft_w_360grid_fortcode, power_360grid_fortcode]: 
        # executable file 
        fort_exe = 'exe/'+code.rsplit('.')[0]+'.exe'
        # compile command
        compile_cmd = ' '.join(['ifort -O3 -o', ldg_code_dir+fort_exe, ldg_code_dir+code, 
            '-L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw'])
        print compile_cmd
        # call compile command 
        subprocess.call(compile_cmd.split())
    #------------------------------------------------------------------------------------------------------------------
    # Generate test MOCKS --------------------------------------------------------------------------------------------
    # randomly displace a number of galaxies
    file_indices = []
    original_files = []
    displaced_files = []
    for i_mock in range(1, n_mocks+1): 
        for letter in ['a', 'b', 'c', 'd']: 
            file_indices.append([i_mock, letter])
            original_files.append(''.join(['/mount/chichipio2/rs123/MOCKS/LRGFull_zm_geo/gaussian/zspace/', 
                mock_prefix, str(100+i_mock)[1:3], letter, mock_suffix]))
            displaced_files.append(''.join([get_fibcoll_dir('data', catalog='lasdamasgeo'), 
                mock_prefix, str(100+i_mock)[1:3], letter, test_suffix]))

    for i_file, displaced_file in enumerate(displaced_files):  
        if os.path.isfile(displaced_file) == False:  
            displaced_cmd = ' '.join([ldg_code_dir+'exe/'+peakcorr_fortcode.rsplit('.')[0]+'.exe', 
                            nbar_file, original_files[i_file], displaced_file])
            print "Displacing the Mocks!"
            subprocess.call(displaced_cmd.split())
        else: 
            print displaced_file, ' already exists'

    # FFT ---------------------------------------------------------------------------------------------------------
    FFT_exe = ldg_code_dir+'exe/'+fft_w_360grid_fortcode.rsplit('.')[0]+'.exe'

    # FFT for random 
    fft_rand_file = ''.join([get_fibcoll_dir('fft', catalog='lasdamasgeo'), 
        'FFT_', rand_file.rsplit('/')[-1], '.randdisptest.grid', grid, '.P0', str(P0), '.box', box])
    FFT_rand_cmd = ' '.join([FFT_exe, str(Rbox), "1", str(P0), rand_file, fft_rand_file])
    if os.path.isfile(fft_rand_file) == False: 
        print "Building ", fft_rand_file 
        subprocess.call(FFT_rand_cmd.split())
    else: 
        print fft_rand_file, " already exists" 
    
    # FFT for mocks
    fft_files = [] 
    for i_file, file_index in enumerate(file_indices): 
        fft_file = ''.join([get_fibcoll_dir('fft', catalog='lasdamasgeo'), 
            'FFT_', displaced_files[i_file].rsplit('/')[-1], '.grid', grid, '.P0', str(P0), '.box', box])
        fft_files.append(fft_file) 
        
        # FFT command
        print 'Building ', fft_file 
        FFT_cmd = ' '.join([FFT_exe, str(Rbox), "0", str(P0), displaced_files[i_file], fft_file]) 
        subprocess.call(FFT_cmd.split()) 

    # P(k) ---------------------------------------------------------------------------------------------------------
    power_exe = ldg_code_dir+'exe/'+power_360grid_fortcode.rsplit('.')[0]+'.exe'
    for i_file, file_index in enumerate(file_indices): 
        power_file = ''.join([get_fibcoll_dir('power', catalog='lasdamasgeo'), 
            'power_', displaced_files[i_file].rsplit('/')[-1], '.randdisttest.grid', grid, '.P0', str(P0), '.box', box])
        
        power_cmd = ' '.join([power_exe, fft_files[i_file], fft_rand_file, power_file, str(sscale)]) 
        print 'Building ', power_file 
        subprocess.call(power_cmd.split())

def build_ldg_peak_corrected_randoms(sigma, fpeak, corr='peak', sanitycheck=False): 
    # Build Corrected Randoms for Las Damas Geo 
    build_corrected_randoms(catalog='LasDamasGeo', 
            catalog_param_list={'n_mock':[[i,j] for i in range(1,41) for j in ['a', 'b', 'c', 'd']]}, 
            corr=corr, corr_param={'sigma':sigma, 'fpeak':fpeak}, sanitycheck=sanitycheck)

def build_ldg_upweight_corrected_randoms(): 
    # Build Corrected Randoms for Las Damas Geo 
    build_corrected_randoms(catalog='LasDamasGeo', 
            catalog_param_list={'n_mock':[[i,j] for i in range(1,41) for j in ['a', 'b', 'c', 'd']]}, 
            corr='upweight')

def build_ldg_peak_corrected_nbar(): 
    ''' 
    Build peak corrected nbar for Las Damas Geo
    '''
    for fpeak in [0.1*i for i in range(0,11)]: 
        ldg_peak_corr_nbar = nbar(catalog='lasdamasgeo', corr='peaknbar', corr_param={'sigma':5.3, 'fpeak':fpeak}) 
        ldg_peak_corr_nbar.writenbar()

if __name__=="__main__":
    '''
    #lasdamasgeo_fibcoll_pk_rand_disp_test(1)
    ldg_catalog_list = [[i,j] for i in range(1,2) for j in ['a', 'b', 'c', 'd']]
    #build_peakcorr_fibcoll(catalog='LasDamasGeo', catalog_param_list={'n_mock':ldg_catalog_list}, corr_param={'sigma':5.3, 'fpeak':0.6})
    lasdamasgeo_fibcoll_pk(1, corr='peaktest')
    '''
    plot_bispec_fibcolcorr_comparison(BorQ='B', x_axis='triangles', triangle='all', 
        catalog='PTHalo', filespec={'version':'v11p0', 'n_mock':range(1,11), 'Nrandom':100}, 
        correction_method=['true', 'delta'], resid='True')

    plot_bispec_fibcolcorr_comparison(BorQ='B', x_axis='triangles', triangle='acute', 
        catalog='PTHalo', filespec={'version':'v11p0', 'n_mock':range(1,11), 'Nrandom':100}, 
        correction_method=['true', 'delta'], resid='True')

    plot_bispec_fibcolcorr_comparison(BorQ='B', x_axis='triangles', triangle='obtuse', 
        catalog='PTHalo', filespec={'version':'v11p0', 'n_mock':range(1,11), 'Nrandom':100}, 
        correction_method=['true', 'delta'], resid='True')

    plot_bispec_fibcolcorr_comparison(BorQ='B', x_axis='triangles', triangle='extended', 
        catalog='PTHalo', filespec={'version':'v11p0', 'n_mock':range(1,11), 'Nrandom':100}, 
        correction_method=['true', 'delta'], resid='True')
    '''
    corr_method = ['true', 'peak', 'randdisptest'] 
    sigmas = [0, 5.3, 5.3]
    fpeaks = [0, 0.1, 0.5]

    plot_powerspec_fibcolcorr_comparison(catalog='lasdamasgeo', filespec={'n_mock':ldg_catalog_list}, 
            correction_method=corr_method, corrspec={'sigma':sigmas, 'fpeak':fpeaks}, resid='False')
    plot_powerspec_fibcolcorr_comparison(catalog='lasdamasgeo', filespec={'n_mock':ldg_catalog_list}, 
            correction_method=corr_method, corrspec={'sigma':sigmas, 'fpeak':fpeaks}, resid='True')
    #build_ldg_peak_corrected_randoms() 
    '''
