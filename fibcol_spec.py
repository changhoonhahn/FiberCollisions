import fibcol_data as fc_data
import numpy as np
import os.path
import subprocess
import cosmolopy as cosmos
import fibcol_fft as fc_fft
import fibcol_utility as fc_util

# Classes ------------------------------------------------------------
class Spec: 
    def __init__(self, spectrum, Igal_Irand=False, **cat_corr):
        '''
        bispectrum file 
        specify catalog, version, mock file number, file specifications (e.g. Nrandom), 
        fiber collision correction method, correction specifications (e.g. sigma, fpeak)
        '''
        catalog = cat_corr['catalog'] 
        correction = cat_corr['correction'] 
        spec = cat_corr['spec']
        self.cat_corr = cat_corr        # store the catalogue/correction/spec dictionary

        # powerspectrum or bispectrum
        self.spectrum = spectrum
        spec_dir = '/mount/riachuelo1/hahn/'+spectrum.lower()+'/'
        if spectrum.lower() == 'power':                         # set flags
            spec_file_flag = 'power_'
        elif spectrum.lower() == 'bispec': 
            spec_file_flag = 'bisp_'
        
        try: 
            spec['quad']
        except KeyError: 
            spec['quad'] = False

        if spec['quad'] == True:                                # for Quadrupole code 
            spec_file_flag = spec_file_flag+'quad_'

        # Tiling Mock ------------------------------------------------------------------------------------------------
        if catalog['name'].lower() == 'tilingmock':  
            # Tiling Mock directory 
            file_dir = spec_dir+'tiling_mocks/'      

            data = fc_data.galaxy_data('data', readdata=False, **cat_corr)
            data_file = data.file_name
            if correction['name'].lower() == 'shotnoise': data_file = data_file+'.shotnoise'
            if correction['name'].lower() == 'floriansn': data_file = data_file+'.floriansn'
            if correction['name'].lower() == 'hectorsn': data_file = data_file+'.hectorsn'

            # file naem  
            file_prefix = spec_file_flag+data_file.rsplit('/')[-1]+'.corrnbar'

            # file ending  
            if spectrum == 'bispec': 
                file_suffix = '.grid360.nmax.nstep3.P020000.box3600'
            elif spectrum == 'power': 
                file_suffix = '.grid'+str(spec['grid'])+'.P0'+str(spec['P0'])+'.box'+str(spec['box'])
            
            # specify correction 
            file_corr = ''          # specified in file_prefix within the data file name 

            self.scale = spec['box']
        elif catalog['name'].lower() == 'cmass': 
            file_dir = spec_dir

            file_prefix = 'power-cmass-dr12v4-N-Reid-weights-zlim-ngalsys-3600lbox-360grid-180bin.dat'
            file_corr = ''
            file_suffix = ''

        # QPM ------------------------------------------------------------------------------------------------------------
        elif catalog['name'].lower() == 'qpm': 
            file_dir = spec_dir+'QPM/dr12d/'                       # QPM directory 

            data = fc_data.galaxy_data('data', readdata=False, **cat_corr)
            data_file = data.file_name 

            # file naem  
            file_prefix = spec_file_flag+data_file.rsplit('/')[-1]
            
            if correction['name'].lower() == 'shotnoise': file_prefix = file_prefix+'.shotnoise'
            if correction['name'].lower() == 'floriansn': file_prefix = file_prefix+'.floriansn'
            if correction['name'].lower() == 'hectorsn': file_prefix = file_prefix+'.hectorsn'

            # file ending  
            if spectrum == 'bispec': 
                file_suffix = '.grid360.nmax.nstep3.P020000.box3600'
            elif spectrum == 'power': 
                file_suffix = '.grid'+str(spec['grid'])+'.P0'+str(spec['P0'])+'.box'+str(spec['box'])
            
            # specify correction 
            file_corr = ''          # specified in file_prefix within the data file name 

            self.scale = spec['box']

        # LasDamas Geomtry Mocks ----------------------------------------------------------------------------------
        elif catalog['name'].lower() == "lasdamasgeo": 
            # Las Damas Directory
            file_dir = ''.join([spec_dir, '/LasDamas/Geo/'])

            # e.g. power_sdssmock_gamma_lrgFull_zm_oriana19a_no.rdcz.dat.grid360.P020000.box3600
            # file beginning
            data = fc_data.galaxy_data('data', readdata=False, **cat_corr)
            data_file = data.file_name

            file_prefix = spec_file_flag+data_file.rsplit('/')[-1]

            if correction['name'].lower() == 'shotnoise': file_prefix = file_prefix+'.shotnoise'
            if correction['name'].lower() == 'floriansn': file_prefix = file_prefix+'.floriansn'
            if correction['name'].lower() == 'hectorsn': file_prefix = file_prefix+'.hectorsn'

            # file ending  
            if spectrum == 'bispec': 
                file_suffix = '.grid360.nmax.nstep3.P020000.box3600'
            elif spectrum == 'power': 
                file_suffix = '.grid'+str(spec['grid'])+'.P0'+str(spec['P0'])+'.box'+str(spec['box'])
            
            # specify correction 
            file_corr = ''          # specified in file_prefix within the data file name 

            # survey scale 
            self.scale = spec['box']
        
            ''' 
            # correction specifier 
            if correction['name'].lower() in ('true'): 
                # no correction
                file_corr = 'dat'
            elif correction['name'].lower() in ('shotnoise', 'hectorsn', 'floriansn'):
                file_corr = ''
            # upweight correction -------------------------------------------------------------------------------------
            elif correction['name'].lower() in ('delta', 'upweight'): 
                # correction flag
                file_corr = 'fibcoll.dat.upweight'

            # correction with peak and tail ---------------------------------------------------------------------------
            elif correction['name'].lower() in ('peak', 'peaknbar', 'peaktest', 'peakshot', 'vlospeakshot') : 
                
                if correction['name'].lower() == 'peak':
                    # Correction for poor naming convention 
                    correction['name'] = 'peaknbar'

                # specify peak correction fit (expon or gauss) 
                if correction['fit'].lower() in ('gauss', 'expon'): 
                    fit_str = correction['fit'].lower()

                    # correction flag 
                    file_corr = 'fibcoll.dat.'+fit_str+'.'+correction['name'].lower()+\
                            '.sigma'+str(correction['sigma'])+'.fpeak'+str(correction['fpeak'])+'.corrnbar'

                elif correction['fit'].lower() in ('real'): 
                    fit_str = correction['fit'].lower()
                    # correction flag 
                    file_corr = 'fibcoll.dat.'+fit_str+'.'+correction['name'].lower()+'.fpeak'+str(correction['fpeak'])+'.corrnbar'

                else: 
                    raise NameError('peak fit has to be specified: gauss or expon') 


            # correction with *only* peak --------------------------------------------------------------------------------
            elif correction['name'].lower() in ('allpeak', 'allpeakshot'):

                # specify peak correction fit (expon or gauss) 
                if correction['fit'].lower() == 'gauss': 
                    fit_str = 'gauss' 
                elif correction['fit'].lower() == 'expon': 
                    fit_str = 'expon'
                else: 
                    raise NameError('peak fit has to be specified: gauss or expon') 

                # correction flag 
                file_corr = 'fibcoll.dat.'+fit_str+'.'+correction['name'].lower()+\
                        '.sigma'+str(correction['sigma'])+'.fpeak'+str(correction['fpeak'])+'.corrnbar'
            
            # Test adjustments --------------------------------------------------------------------------------------------
            elif (correction['name'].lower() == 'randrm'): 

                # correction flag 
                file_corr = '.dat.'+correction['name'].lower()+'.corrnbar'

            else: 
                raise NameError('Correction method not supported') 
            '''

        # combine to make file
        self.file_name = ''.join([file_dir, file_prefix, file_corr, file_suffix]) # combine file parts 
        #print self.file_name

    def readfile(self): 
        '''
        Read power/bi-spectrum file and import values of interest  
        '''
        correction = self.cat_corr['correction']
        spec = self.cat_corr['spec']

        # read file 
        file = np.loadtxt(self.file_name) 
        
        # Powerspectrum ----------------------------------------------------------------------------------------------------
        if self.spectrum == 'power':
            # columns of data to be imported
            self.columns = ['k', 'Pk']
            self.k  = file[:,0]
            if spec['quad'] == True: 
                if correction['name'].lower() in ('peakshot', 'vlospeakshot'): 
                    self.Pk = file[:,5]
                else: 
                    self.Pk = file[:,1]
            else: 
                self.Pk = file[:,1]
            
            if spec['quad'] == True: 
                self.P2k = file[:,2]            # power spectrum quadrupole
            else: 
                pass

        # Bispectrum ----------------------------------------------------------------------------------------------------
        elif self.spectrum == 'bispec': 
            # columns of data to be imported
            self.columns = ['kfund', 'i_triangles', 'k1', 'k2', 'k3', 'Pk1', 'Pk2', 'Pk3', 
                    'Bk', 'Q', 'avgk', 'kmax']

            k_fund = (2.0*m.pi)/np.float(self.scale)
            self.kfund = k_fund     # fundamental k 

            # import following values from file: 
            self.i_triangle = range(len(file[:,0]))     # triangle numbers 

            # triangle k values  
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

def build_fibcol_power(**cat_corr): 
    '''
    Given catalog_correction dictionary, construct FFT file 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    spec = cat_corr['spec'] 
    
    # data file name 
    data = fc_data.galaxy_data('data', readdata=False, **cat_corr) 

    # FFt file names 
    fft_file = fc_fft.get_fibcol_fft_file('data', **cat_corr) 
    fft_rand_file = fc_fft.get_fibcol_fft_file('random', **cat_corr) 
    
    power = Spec('power', **cat_corr) 
    power_file = power.file_name 
    
    if spec['quad'] == True:                # for quadrupole power
        power_code = fc_util.fortran_code('quadpower', **cat_corr) 
    else: 
        power_code = fc_util.fortran_code('power', **cat_corr) 
    power_exe = fc_util.fortran_code2exe(power_code)
    
    # code and exe modification time 
    power_code_mod_time = os.path.getmtime(power_code)
    if os.path.isfile(power_exe) == False: 
        power_exe_mod_time = 0 
    else: 
        power_exe_mod_time = os.path.getmtime(power_exe)

    # if code was changed since exe file was last compiled then 
    # compile power code 
    if power_exe_mod_time < power_code_mod_time: 
        fc_util.compile_fortran_code(power_code) 
    
    # LasDamas Geo and Tiling Mock ----------------------------------------------------------------
    if catalog['name'].lower() in ('lasdamasgeo', 'tilingmock', 'qpm'): 
        if spec['quad'] == True:            # for quadrupole code NOTE: ORDER OF RAND AND MOCK FILES ARE REVERSED
            power_cmd = ' '.join([power_exe, fft_rand_file, fft_file, power_file, 
                str(spec['sscale']), str(spec['grid']/2)]) 
        else: 
            power_cmd = ' '.join([power_exe, fft_file, fft_rand_file, power_file, str(spec['sscale'])]) 
        print power_cmd
        subprocess.call(power_cmd.split()) 
    else: 
        raise NameError('not yet coded') 
            
    return power_file  

"""         #COMMENTED AWAY, UNRELIABLE 
def build_fibcol_power_Igal_Iran(**cat_corr):           
    '''
    Given catalog_correction dictionary, construct FFT file 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    spec = cat_corr['spec'] 
    
    # data file name 
    data = fc_data.galaxy_data('data', readdata=False, **cat_corr) 

    # FFt file names 
    fft_file = get_fibcol_fft_file('data', **cat_corr) 
    fft_rand_file = get_fibcol_fft_file('random', **cat_corr) 

    power = fc_spec.spec('power', **cat_corr) 
    power_file = power.file_name+'_Igal_Iran'

    power_code = '/home/users/hahn/powercode/FiberCollision/LasDamas/Geo/power_ldg_fkp_Igal_Iran_360grid_180bin.f'
    power_exe = fc_util.fortran_code2exe(power_code)
   
    '''
    # code and exe modification time 
    power_code_mod_time = time.ctime(os.path.getmtime(power_code))
    power_exe_mod_time = time.ctime(os.path.getmtime(power_exe))

    # if code was changed since exe file was last compiled then 
    # compile power code 
    if (time.ctime(os.path.getmtime(power_exe))< power_code_mod_time) or (os.path.isfile(power_exe) == False): 
    '''
    fc_util.compile_fortran_code(power_code) 

    if catalog['name'].lower() == 'lasdamasgeo': 
        power_cmd = ' '.join([power_exe, fft_file, fft_rand_file, power_file, str(spec['sscale'])]) 
        print power_cmd
        subprocess.call(power_cmd.split()) 
            
    return power_file  
"""

def build_fibcol_bispec(**cat_corr): 
    '''
    Given catalog_correction dictionary, construct bispec file 
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    spec = cat_corr['spec'] 
    
    # data file name 
    data = fc_data.galaxy_data('data', readdata=False, **cat_corr) 

    # FFt file names 
    fft_file = get_fibcol_fft_file('data', **cat_corr) 
    fft_rand_file = get_fibcol_fft_file('random', **cat_corr) 

    bispec = fc_spec.spec('bispec', **cat_corr) 
    bispec_file = bispec.file_name 

    bispec_code = fc_util.fortran_code('bispec', **cat_corr) 
    bispec_exe = fc_util.fortran_code2exe(power_code)
    
    # code and exe modification time 
    power_code_mod_time = time.ctime(os.path.getmtime(power_code))
    power_exe_mod_time = time.ctime(os.path.getmtime(power_exe))

    # if code was changed since exe file was last compiled then 
    # compile power code 
    if (power_exe_mod_time < power_code_mod_time) or (os.path.isfile(bispec_exe) == False): 
        fc_util.compile_fortran_code(bispec_code) 

    if catalog['name'].lower() == 'lasdamasgeo': 
        bispec_cmd = ' '.join([bispec_exe, '2', fft_rand_file, fft_file, bispec_file]) 
        print power_cmd
        subprocess.call(power_cmd.split()) 
            
    return power_file  
