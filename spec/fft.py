import numpy as np 
import os.path
import time 
import subprocess
import cosmolopy as cosmos

from data import Data 
from util.direc import direc
from util.fortran import Fcode

class Fft(object): 

    def __init__(self, DorR, cat_corr, **kwargs): 
        """ A class that describes the FFT of galaxy simulated/observed data 
        """
        if 'spec' not in cat_corr.keys(): 
            # default spectrum parameters
            cat_corr['spec'] = {
                    'P0': 20000, #P0 
                    'Lbox': 3600.0, 
                    'Ngrid':360, 
                    'quad': False
                    }

        self.cat_corr = cat_corr 
        self.kwargs = kwargs
        self.type = DorR

        self.file_name = self.file()

    def file(self): 
        """ FFT data file name 
        """

        corrdict = (self.cat_corr)['correction']
        specdict = (self.cat_corr)['spec'] 
    
        fft_dir = direc('fft', self.cat_corr)

        galdata = Data(self.type, self.cat_corr, **self.kwargs)  # data class 
        self.data_file = galdata.file_name # galaxy data file

        if not os.path.isfile(self.data_file): 
            raise ValueError()
    
        # FFT label 
        try: 
            if specdict['quad']: 
                fft_str = 'FFT_Q_'
            else: 
                fft_str = 'FFT_'
        except KeyError: 
            fft_str = 'FFT_'
    
        fft_corr_str = ''
        if (corrdict['name'].lower() in ('floriansn', 'hectorsn')) & (self.type != 'random'):
            fft_corr_str = ''.join(['.', corrdict['name'].lower()])

        # FFTs from data file 
        fft_file = ''.join([
            fft_dir, 
            fft_str, (self.data_file).rsplit('/')[-1], 
            fft_corr_str,
            '.grid', str(specdict['Ngrid']), 
            '.P0', str(specdict['P0']), 
            '.box', str(specdict['Lbox'])
            ])

        return fft_file  

    def build(self): 
        """ Run FFT FORTRAN code to calculate FFT of data
        """
        
        catdict = (self.cat_corr)['catalog']
        corrdict = (self.cat_corr)['correction']
        specdict = (self.cat_corr)['spec'] 
        
        if not os.path.isfile(self.data_file):
            galdata = Data(self.type, self.cat_corr, **self.kwargs) 
            galdata.build()

        if 'quad' not in specdict.keys(): 
            specdict['quad'] = False

        if not specdict['quad']:       # quadrupole or regular FFT code
            fft_type = 'fft'
        else:  
            NotImplementedError()

        codeclass = Fcode(fft_type, self.cat_corr) 
        fftcode = codeclass.code
        fftexe = codeclass.fexe()
        
        # code and exe modification time 
        fftcode_t_mod, fftexe_t_mod = codeclass.mod_time()

        if fftexe_t_mod < fftcode_t_mod: 
            codeclass.compile()
                
        fft_file = self.file() 
        if self.type == 'data': 
            n_DorR = 0
        elif self.type == 'random': 
            n_DorR = 1

        if not specdict['quad']:       # NOT Quadrupole
            
            # command line call 
            FFTcmd = codeclass.commandline_call(
                    DorR = self.type, 
                    datafile = self.data_file,
                    fftfile = self.file_name
                    ) 

            if 'clobber' not in (self.kwargs).keys(): 
                bool_clobber = False

            if any([not os.path.isfile(self.file_name), bool_clobber]):
                print ''
                print '-----------------------'
                print 'Constructing '
                print self.file_name  
                print '-----------------------'
                print ''
                print FFTcmd
                print '-----------------------'

                subprocess.call(FFTcmd.split())
            else: 
                print ''
                print '-----------------------'
                print self.file_name  
                print 'Already Exists'
                print '-----------------------'
                print ''

        else:       
            raise NotImplementedError()
        return None 

if __name__=='__main__': 
    cat_corr = {'catalog': {'name': 'nseries', 'n_mock': 1}, 'correction': {'name': 'upweight'}}
    for DorR in ['data', 'random']:
        fftee = Fft(DorR, cat_corr, clobber=True)
        print fftee.build()

"""
        # Quad FFT argument sequence (SUBJECT TO CHANGE) 

            # determine "idata"
            if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz'): 
                idata = 2
                ifc = 0 
            elif catalog['name'].lower() == 'qpm': 
                idata = 3 
                ifc = 0 
            elif catalog['name'].lower() == 'tilingmock':
                idata = 9 
                ifc = 0 
            elif catalog['name'].lower() == 'nseries': 
                idata = 10 
                ifc = 0 
            elif catalog['name'].lower() == 'ldgdownnz':  
                idata = 11 
                ifc = 0 
            elif 'bigmd' in catalog['name'].lower(): 
                idata = 12
                ifc = 0 
            elif catalog['name'].lower() == 'cmass': 
                idata = 13 
                ifc = 0 
            else: 
                raise NameError('not included in Quadrupole code') 
                
            # bash commend 
            # is of the form 
            # FFT_FKP_BOSS_cic_il4_v3.exe idata box Ngrid interpol iflag P0  ifc icomp input_file output_file
            # icomp is hardcoded 0 so that it takes into account completeness!
            FFT_cmd = ' '.join([
                FFT_exe, str(idata), 
                str(spec['box']), str(spec['grid']), 
                "4", str(DorR_number), str(spec['P0']), 
                str(ifc), "0", data_file, fft_file]) 
            print FFT_cmd

            if DorR.lower() == 'data':  # don't bother checking if the file exists for mocks and run the damn thing 
                subprocess.call(FFT_cmd.split()) 

            elif DorR.lower() == 'random':      # random takes longer so check to see if it exists first
                # call FFT randomc ommand 
                if os.path.isfile(fft_file) == False: 
                    print "Building ", fft_file 
                    subprocess.call(FFT_cmd.split())
                else: 
                    print fft_file, " already exists" 

            print 'Constructing ', 
"""
