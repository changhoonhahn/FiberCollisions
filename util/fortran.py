'''

Interface with FORTRAN code 

'''
import subprocess
import os.path

class Fcode: 

    def __init__(self, type, cat_corr): 
        """ Class to describe FORTRAN code for fiber collision corrected 
        powerspectrum calculations 

        """
        self.cat_corr = cat_corr
        self.type = type
    
        fcode_dir = '/home/users/hahn/powercode/FiberCollisions/fortran/'
        
        if type == 'fft':       # fft code
            f_name = 'FFT_fkp.f'
        elif type == 'pk':      # P(k) code
            f_name = 'power_IgalIrand.f' 
        else: 
            raise NotImplementedError()
        
        f_code = ''.join([fcode_dir, f_name])
        self.code = f_code
        self.exe = None 

    def fexe(self): 
        """ Fortran executable that corresponds to fortran code
        """
        code_dir = ''.join(['/'.join((self.code).split('/')[:-1]), '/'])
        code_file = (self.code).split('/')[-1]
    
        fort_exe = ''.join([code_dir, 'exe/', '.'.join(code_file.rsplit('.')[:-1]), '.exe'])
        self.exe = fort_exe 
        
        return fort_exe 

    def compile(self):
        """ Compile fortran code
        """

        fort_exe = self.fexe() 

        # compile command
        if self.code == '/home/users/hahn/powercode/FiberCollisions/fortran/FFT_FKP_BOSS_cic_il4_v3.f': 
            compile_cmd = ' '.join(['ifort -fast -o', fort_exe, self.code, '-L/usr/local/fftw_intel_s/lib -lsrfftw -lsfftw -lm'])
        elif self.code == '/home/users/hahn/powercode/FiberCollisions/fortran/power_FKP_SDSS_BOSS_v3.f': 
            compile_cmd = ' '.join(['ifort -fast -o', fort_exe, self.code])
        else: 
            compile_cmd = ' '.join(['ifort -O3 -o', fort_exe, self.code, '-L/usr/local/fftw_intel_s/lib -lsfftw -lsfftw'])
        print ' ' 
        print 'Compiling -----'
        print compile_cmd
        print '----------------'
        print ' ' 

        # call compile command 
        subprocess.call(compile_cmd.split())

        return None 

    def commandline_call(self, **kwargs): 
        """ Command line call for Fortran code
        """
        
        fcode_t_mod, fexe_t_mod = self.mod_time()
        if fcode_t_mod < fcode_t_mod: 
            raise ValueError("Compile failed")

        fort_exe = self.fexe() 

        if self.type == 'fft': 

            if not all([kwarg in kwargs.keys() for kwarg in ['DorR', 'datafile', 'fftfile']]):
                err_msg = ''.join(['For type = ', self.type, " 'DorR', 'datafile', 'fftfile' must be specified in kwargs"])
                raise KeyError(err_msg)
            
            catname = ((self.cat_corr)['catalog'])['name']
            specdict = (self.cat_corr)['spec']

            if catname == 'nseries': 
                n_cat = 7
            else: 
                raise NotImplementedError()
            
            try: 
                if kwargs['cosmology'] == 'survey': 
                    n_cosmo = 1
                else: 
                    n_cosmo = 0 
            except KeyError: 
                n_cosmo = 0 

            if kwargs['DorR'] == 'data': 
                n_DorR = 0 
            elif kwargs['DorR'] == 'random': 
                n_DorR = 1 

            cmdline_call = ' '.join([
                self.exe, 
                str(n_cat), 
                str(n_cosmo), 
                str(specdict['box']), 
                str(specdict['grid']), 
                str(n_DorR),
                str(specdict['P0']), 
                kwargs['datafile'], 
                kwargs['fftfile']
                ])
        else: 
            raise NotImplementError()

        return cmdline_call

    def mod_time(self): 
        """ Modification time of .f and .exe file 
        """
        
        if self.exe == None: 
            self.fexe()
        
        if not os.path.isfile(self.code): 
            fcode_t_mod = 0 
        else: 
            fcode_t_mod = os.path.getmtime(self.code)

        if not os.path.isfile(self.exe): 
            fexe_t_mod = 0 
        else: 
            fexe_t_mod = os.path.getmtime(self.exe)
        
        return [fcode_t_mod, fexe_t_mod]

if __name__=="__main__": 
    cat_corr = {'catalog': {'name': 'nseries', 'n_mock': 1}, 'correction': {'name': 'upweight'}}
    code = fcode('fft', cat_corr)
    print code.code
    print code.fexe()
    print code.compile()

"""
        if fft_power.lower() == 'fft':          # FFT code
            if correction['name'].lower() == 'floriansn': 
                code_file = 'FFT-nseries-fkp-w-nbar-florian-'+str(spec['grid'])+'grid.f'
            elif correction['name'].lower() == 'hectorsn': 
                code_file = 'FFT-nseries-fkp-w-nbar-hector-'+str(spec['grid'])+'grid.f'
            else: 
                code_file = 'FFT-nseries-fkp-w-nbar-'+str(spec['grid'])+'grid.f'
    
        elif fft_power.lower() == 'power':      # power code
            if correction['name'].lower() in ('true', 'upweight', 'peaknbar'):
                # normal FKP shot noise correction
                if spec['grid'] == 360: 
                    code_file = 'power-nseries-fkp-w-nbar-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    code_file = 'power-nseries-fkp-w-nbar-960grid-480bin.f'
            elif correction['name'].lower() in \
                    ('peakshot', 'photozpeakshot', 'photozenvpeakshot', 
                            'shotnoise', 'floriansn', 'noweight', 'hectorsn', 'peakshot_dnn'): 
                # FKP with Igal Irand shot noise correction 
                if spec['grid'] == 360: 
                    code_file = 'power-nseries-fkp-w-nbar-Igal-Irand-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    code_file = 'power-nseries-fkp-w-nbar-Igal-Irand-960grid-480bin.f'
            elif 'scratch' in correction['name'].lower(): 
                # FKP with Igal Irand shot noise correction (for scratch pad corrections) 
                if spec['grid'] == 360: 
                    code_file = 'power-nseries-fkp-w-nbar-Igal-Irand-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    code_file = 'power-nseries-fkp-w-nbar-Igal-Irand-960grid-480bin.f'
            else: 
                raise NameError('asldkfjasdf') 

        # quadrupole codes --------------------------------------------
        # regardess of catalog or correction TEMPORARILY HARDCODED HERE FOR TEST RUN 
        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            code_file = 'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            code_file = 'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError("not Yet coded") 
        
        f_code = ''.join([code_dir, code_file]) 

    if catalog['name'].lower() == 'lasdamasgeo':            # LasDamasGeo ----------------------

        cat_dir = 'LasDamas/Geo/' 
       
        if fft_power.lower() == 'fft':                      # FFT -----------------
    
            if correction['name'].lower() == 'floriansn':     # Beutler+2014
                f_code = ''.join(['FFT_ldg_fkp_w_florian_', str(spec['grid']), 'grid.f']) 
            elif correction['name'].lower() == 'hectorsn':  # Gil-Marin+2014
                f_code = ''.join(['FFT_ldg_fkp_w_hector_', str(spec['grid']), 'grid.f'])
            else: 
                f_code = ''.join(['FFT_ldg_fkp_w_', str(spec['grid']), 'grid.f'])

        elif fft_power.lower() == 'power':                  # power ----------------

            if correction['name'].lower() in ('true', 'upweight', 'peaknbar', 'bigfc'):  
                # FKP estimator
                if spec['grid'] == 360: 
                    f_code = ldg_code_dir+'power_ldg_fkp_360grid_180bin.f'
                elif spec['grid'] == 960: 
                    f_code = ldg_code_dir+'power_ldg_fkp_960grid_480bin.f'

            elif correction['name'].lower() in ('peakshot', 'allpeakshot', 'noweight', 
                    'shotnoise', 'vlospeakshot', 'floriansn', 'hectorsn', 'peakshot_dnn', 
                    'bigfc_peakshot'): 
                # Igal+Irand shot noise incorporated
                if spec['grid'] == 360: 
                    f_code = ldg_code_dir+'power_ldg_fkp_Igal_Iran_360grid_180bin.f'
                elif spec['grid'] == 960: 
                    f_code = ldg_code_dir+'power_ldg_fkp_Igal_Iran_960grid_480bin.f'
            else: 
                raise NameError('what?')
    
        # quadrupole codes ------------------------------------------------
        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError('asdflkajsdf') 
    
    elif catalog['name'].lower() == 'ldgdownnz':            # LasDamasGeo downsampled ---------
        ldg_code_dir = '/home/users/hahn/powercode/FiberCollisions/LasDamas/Geo/' 
       
        if fft_power.lower() == 'fft':      # FFT -----------------
            # only one implemented so far 
            f_code = ldg_code_dir+'FFT_ldg_fkp_w_down_nz_'+str(spec['grid'])+'grid.f'

        elif fft_power.lower() == 'power':  # power ----------------

            if correction['name'].lower() in ('true', 'upweight', 'peaknbar', 'bigfc'):  
                # Original FKP estimator
                if spec['grid'] == 360: 
                    f_code = ldg_code_dir+'power_ldg_fkp_360grid_180bin.f'
                elif spec['grid'] == 960: 
                    f_code = ldg_code_dir+'power_ldg_fkp_960grid_480bin.f'

            elif correction['name'].lower() in ('peakshot', 'bigfc_peakshot'): 
                # Igal+Irand shot noise incorporated
                if spec['grid'] == 360: 
                    f_code = ldg_code_dir+'power_ldg_fkp_Igal_Iran_360grid_180bin.f'
                elif spec['grid'] == 960: 
                    f_code = ldg_code_dir+'power_ldg_fkp_Igal_Iran_960grid_480bin.f'

            else: 
                raise NotImplementedError('what?')

        elif fft_power.lower() == 'quadfft':  # Quadrupole FFT
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 

        elif fft_power.lower() == 'quadpower': # Quadrupole power
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'

        else: 
            raise NameError('asdflkajsdf') 
    
    elif catalog['name'].lower() == 'tilingmock':               # Tiling Mock ----------------
        code_dir = '/home/users/hahn/powercode/FiberCollisions/TilingMock/'
    
        if fft_power.lower() == 'fft':
            if correction['name'].lower() == 'floriansn': 
                f_code = code_dir+'FFT-fkp-tm-w-nbar-florian-'+str(spec['grid'])+'grid.f'
            elif correction['name'].lower() == 'hectorsn': 
                f_code = code_dir+'FFT-fkp-tm-w-nbar-hector-'+str(spec['grid'])+'grid.f'
            else: 
                f_code = code_dir+'FFT-fkp-tm-w-nbar-'+str(spec['grid'])+'grid.f'

        elif fft_power.lower() == 'power': 
            if correction['name'].lower() in ('peakshot', 'allpeakshot', 'shotnoise', 'floriansn', 'hectorsn', 'vlospeakshot', 'peakshot_dnn'):
                if spec['grid'] == 360: 
                    f_code = code_dir+'power-fkp-tm-w-nbar-Igal-Irand-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    f_code = code_dir+'power-fkp-tm-w-nbar-Igal-Irand-960grid-480bin.f'
            else: 
                if spec['grid'] ==360: 
                    f_code = code_dir+'power-fkp-tm-w-nbar-360grid-180bin.f'
                elif spec['grid'] == 960:
                    f_code = code_dir+'power-fkp-tm-w-nbar-960grid-480bin.f'

        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError('asdlkfjalksdjfklasjf')

    elif catalog['name'].lower() == 'qpm':                  # QPM -----------------------
        code_dir = '/home/users/hahn/powercode/FiberCollisions/QPM/dr12d/'
        
        if fft_power.lower() == 'fft': 
            # FFT code
            if correction['name'].lower() == 'floriansn': 
                f_code = code_dir+'FFT-qpm-fkp-w-nbar-florian-'+str(spec['grid'])+'grid.f'
            elif correction['name'].lower() == 'hectorsn': 
                f_code = code_dir+'FFT-qpm-fkp-w-nbar-hector-'+str(spec['grid'])+'grid.f'
            else: 
                f_code = code_dir+'FFT-qpm-fkp-w-nbar-'+str(spec['grid'])+'grid.f'
    
        elif fft_power.lower() == 'power': 
            if correction['name'].lower() in ('true', 'upweight', 'peaknbar'):
                # normal FKP shot noise correction
                if spec['grid'] == 360: 
                    f_code = code_dir+'power-qpm-fkp-w-nbar-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    f_code = code_dir+'power-qpm-fkp-w-nbar-960grid-480bin.f'

            elif correction['name'].lower() in \
                    ('peakshot', 'shotnoise', 'floriansn', 'noweight', 
                            'hectorsn', 'vlospeakshot', 'peakshot_dnn'): 
                if spec['grid'] == 360: 
                    # Igal Irand shot noise correction 
                    f_code = code_dir+'power-qpm-fkp-w-nbar-Igal-Irand-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    # Igal Irand shot noise correction 
                    f_code = code_dir+'power-qpm-fkp-w-nbar-Igal-Irand-960grid-480bin.f'

        # quadrupole codes --------------------------------------------
        # regardess of catalog or correction TEMPORARILY HARDCODED HERE FOR TEST RUN 
        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError("not Yet coded") 
    
    elif catalog['name'].lower() == 'nseries':              # N series ------------------
        code_dir = 'Nseries/'   # directory
        
        if fft_power.lower() == 'fft':          # FFT code
            if correction['name'].lower() == 'floriansn': 
                code_file = 'FFT-nseries-fkp-w-nbar-florian-'+str(spec['grid'])+'grid.f'
            elif correction['name'].lower() == 'hectorsn': 
                code_file = 'FFT-nseries-fkp-w-nbar-hector-'+str(spec['grid'])+'grid.f'
            else: 
                code_file = 'FFT-nseries-fkp-w-nbar-'+str(spec['grid'])+'grid.f'
    
        elif fft_power.lower() == 'power':      # power code
            if correction['name'].lower() in ('true', 'upweight', 'peaknbar'):
                # normal FKP shot noise correction
                if spec['grid'] == 360: 
                    code_file = 'power-nseries-fkp-w-nbar-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    code_file = 'power-nseries-fkp-w-nbar-960grid-480bin.f'
            elif correction['name'].lower() in \
                    ('peakshot', 'photozpeakshot', 'photozenvpeakshot', 
                            'shotnoise', 'floriansn', 'noweight', 'hectorsn', 'peakshot_dnn'): 
                # FKP with Igal Irand shot noise correction 
                if spec['grid'] == 360: 
                    code_file = 'power-nseries-fkp-w-nbar-Igal-Irand-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    code_file = 'power-nseries-fkp-w-nbar-Igal-Irand-960grid-480bin.f'
            elif 'scratch' in correction['name'].lower(): 
                # FKP with Igal Irand shot noise correction (for scratch pad corrections) 
                if spec['grid'] == 360: 
                    code_file = 'power-nseries-fkp-w-nbar-Igal-Irand-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    code_file = 'power-nseries-fkp-w-nbar-Igal-Irand-960grid-480bin.f'
            else: 
                raise NameError('asldkfjasdf') 

        # quadrupole codes --------------------------------------------
        # regardess of catalog or correction TEMPORARILY HARDCODED HERE FOR TEST RUN 
        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            code_file = 'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            code_file = 'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError("not Yet coded") 
        
        f_code = ''.join([code_dir, code_file]) 

    elif catalog['name'].lower() == 'patchy':               # PATCHY --------------------
        code_dir = '/home/users/hahn/powercode/FiberCollisions/PATCHY/dr12/v6c/'
        
        if fft_power.lower() == 'fft': 
            # FFT code
            if correction['name'].lower() == 'floriansn': 
                f_code = code_dir+'FFT-patchy-fkp-w-nbar-florian-'+str(spec['grid'])+'grid.f'
            elif correction['name'].lower() == 'hectorsn': 
                f_code = code_dir+'FFT-patchy-fkp-w-nbar-hector-'+str(spec['grid'])+'grid.f'
            else: 
                f_code = code_dir+'FFT-patchy-fkp-w-nbar-'+str(spec['grid'])+'grid.f'
    
        elif fft_power.lower() == 'power': 
            if correction['name'].lower() in ('true', 'upweight', 'peaknbar'):
                # normal FKP shot noise correction
                if spec['grid'] == 360: 
                    f_code = code_dir+'power-patchy-fkp-w-nbar-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    pass
                    f_code = code_dir+'power-patchy-fkp-w-nbar-960grid-480bin.f'
            else:
                raise NotImplementedError('not yet implemented')

        # quadrupole codes --------------------------------------------
        # regardess of catalog or correction TEMPORARILY HARDCODED HERE FOR TEST RUN 
        '''
        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError("not Yet coded") 
        '''
    
    elif 'bigmd' in catalog['name'].lower():                # Big MD ------------------------
        code_dir = '/home/users/hahn/powercode/FiberCollisions/BigMD/'
        
        if fft_power.lower() == 'fft':  # FFT
            if correction['name'].lower() == 'floriansn': 
                raise NotImplementedError('asdfklj')
            elif correction['name'].lower() == 'hectorsn': 
                raise NotImplementedError('asdfklj')
            else: 
                f_code = code_dir+'FFT-bigmd-fkp-w-nbar-'+str(spec['grid'])+'grid.f'
    
        elif fft_power.lower() == 'power': 
            if correction['name'].lower() in ('true', 'upweight'):
                # normal FKP shot noise correction
                if spec['grid'] == 360: 
                    f_code = code_dir+'power-bigmd-fkp-w-nbar-360grid-180bin.f'
                elif spec['grid'] == 960: 
                    f_code = code_dir+'power-bigmd-fkp-w-nbar-960grid-480bin.f'
                elif spec['grid'] == 1920: 
                    f_code = code_dir+'power-bigmd-fkp-w-nbar-1920grid-960bin.f'
                else: 
                    raise NotImplementedError('asdfklj')
            else: 
                raise NotImplementedError('asdfklj')

        # quadrupole codes --------------------------------------------
        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 
        elif fft_power.lower() == 'quadpower': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'
        else: 
            raise NameError("not Yet coded") 
    
    elif catalog['name'].lower() == 'cmass':                # CMASS -------------------------
        code_dir = '/home/users/hahn/powercode/FiberCollisions/CMASS/'
        
        if fft_power.lower() == 'fft':  # FFT
            if correction['name'].lower() == 'floriansn': 
                raise NotImplementedError('asdfklj')
            elif correction['name'].lower() == 'hectorsn': 
                raise NotImplementedError('asdfklj')
            else: 
                if catalog['cosmology'].lower() == 'fiducial': 
                    f_code = ''.join([code_dir,
                        'FFT-cmass-', str(spec['grid']), 'grid.f']) 
                else: 
                    raise NotImplementedError('just dont do it') 
    
        elif fft_power.lower() == 'power': 
            if correction['name'].lower() in ('upweight'):  
                # normal FKP shot noise correction
                if spec['grid'] == 360: 
                    f_code = ''.join([code_dir, 
                        'power-cmass-360grid-180bin.f'])
                    #'power_cmass_fkp_Igal_360grid_180bin.f'])

                elif spec['grid'] == 960: 
                    f_code = ''.join([code_dir, 
                        'power-cmass-Igal-960grid-480bin.f'])
                    #    'power-cmass-960grid-480bin.f'])
                
                elif spec['grid'] == 1920: 
                    f_code = ''.join([code_dir, 
                        'power-cmass-1920grid-960bin.f'])

                else: 
                    raise NotImplementedError('asdfklj')

            elif correction['name'].lower() in ('peakshot'): 
                # corrected Igal-alpha*Irand FKP shot noise correction
                if spec['grid'] == 360: 
                    f_code = ''.join([code_dir, 
                        'power-cmass-Igal-360grid-180bin.f'])
                    #'power_cmass_fkp_Igal_360grid_180bin.f'])

                elif spec['grid'] == 960: 
                    f_code = ''.join([code_dir, 
                        'power-cmass-Igal-960grid-480bin.f'])
            else: 
                raise NotImplementedError('asdfklj')

        # quadrupole codes --------------------------------------------
        elif fft_power.lower() == 'quadfft': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'FFT_FKP_BOSS_cic_il4_v3.f' 

        elif fft_power.lower() == 'quadpower': 
            code_dir = '/home/users/hahn/powercode/FiberCollisions/' 
            f_code = code_dir+'power_FKP_SDSS_BOSS_v3.f'

        else: 
            raise NameError("not Yet coded") 
    
    else: 
        raise NaemError('Not coded!') 

"""
