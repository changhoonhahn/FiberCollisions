import numpy as np 
import os.path
import time 
import subprocess
import cosmolopy as cosmos
import fibcol_nbar as fc_nbar
import fibcol_data as fc_data
import fibcol_spec as fc_spec
import fibcol_utility as fc_util
import plot_fibcol as fc_plot

class FFT: 
    def __init__(self, DorR, quad=False, **cat_corr): 
        ''' FFT class to read/write FFT file 
        '''
        catalog = cat_corr['catalog']
        correction = cat_corr['correction']
        spec = cat_corr['spec'] 

        fft_dir = fc_util.get_fibcoll_dir('fft', **cat_corr) # fft directory 

        data_file = fc_data.get_galaxy_data_file(DorR, **cat_corr)  # data file 
    
        if quad == False: 
            FFT_str = 'FFT_'
        else: 
            FFT_str =  'FFT_quad_'
        
        if (correction['name'].lower() in ('floriansn', 'hectorsn')) & \
                (DorR.lower() != 'random'):
            # Florian+ and Hector+ methods require their own random FFTs
            fft_file = ''.join([fft_dir, 
                FFT_str, data_file.rsplit('/')[-1], '.', correction['name'].lower(), 
                '.grid', str(spec['grid']), '.P0', str(spec['P0']), '.box', str(spec['box'])
                ])
        else: 
            # FFTs from data file 
            fft_file = ''.join([fft_dir, 
                FFT_str, data_file.rsplit('/')[-1], 
                '.grid', str(spec['grid']), '.P0', str(spec['P0']), '.box', str(spec['box'])
                ])

        self.file_name = fft_file 

def get_fibcol_fft_file(DorR, **cat_corr): 
    spec = cat_corr['spec']
    fft = FFT(DorR, quad=spec['quad'], **cat_corr) 
    return fft.file_name

def build_fibcol_fft(DorR, **cat_corr): 
    ''' Construct FFT files given catalog_correction dictionary

    Parameters
    ----------
    DorR : 'data'/'random'
    cat_corr : Catalog and Correction dictionary 

    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    spec = cat_corr['spec'] 
    
    # data file name 
    data_file = fc_data.get_galaxy_data_file(DorR, **cat_corr) 
    print data_file

    if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz', 
            'qpm', 'patchy', 'nseries', 'bigmd'):   
        # Las Damas included because it's constant nbar(z)
        pass
    else: 
        data_file = data_file+'.corrnbar'

    try:                # if "quad" is not specified, then it's not used.
        spec['quad'] 
    except KeyError: 
        spec['quad'] = False

    if not spec['quad']:       # quadrupole or regular FFT code
        FFT_code = fc_util.fortran_code('fft', **cat_corr) 
    else:  
        FFT_code = fc_util.fortran_code('quadfft', **cat_corr)  # FFT code 

    FFT_exe = fc_util.fortran_code2exe(FFT_code)            # exe file 
    
    # code and exe modification time 
    FFT_code_mod_time = os.path.getmtime(FFT_code)
    if os.path.isfile(FFT_exe) == False:  
        FFT_exe_mod_time = 0 
    else: 
        FFT_exe_mod_time = os.path.getmtime(FFT_exe)

    # if code was changed since exe file was last compiled then compile fft code 
    if FFT_exe_mod_time < FFT_code_mod_time: 
        fc_util.compile_fortran_code(FFT_code) 
            
    fft_file = get_fibcol_fft_file(DorR, **cat_corr) 
    print 'Constructing ', fft_file 

    if DorR.lower() == 'data': 
        DorR_number = 0
    elif DorR.lower() == 'random': 
        DorR_number = 1

    if not spec['quad']:       # NOT Quadrupole
        if catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz', 
                'tilingmock', 'qpm', 'patchy', 'nseries', 'bigmd'): 
            # Mocks: LasDamas Geo, Tiling Mock, QPM, PATCHY, Nseries, BigMD

            # get bash command 
            FFT_cmd = ' '.join([
                FFT_exe, str(spec['Rbox']), str(DorR_number), str(spec['P0']), 
                data_file, fft_file]) 
            print FFT_cmd
    
            if DorR.lower() == 'data':  
                # don't bother checking if the file exists for mocks and run the damn thing 
                subprocess.call(FFT_cmd.split()) 

            elif DorR.lower() == 'random':      
                # random takes longer so check to see if it exists first
                # call FFT randomc ommand 
                if os.path.isfile(fft_file) == False: 
                    print "Building ", fft_file 
                    subprocess.call(FFT_cmd.split())
                else: 
                    print fft_file, " already exists" 
        else: 
            raise NameError('not yet coded!') 

    else:       # For Quadruopole code ----------------------------------------------------
        # Quad FFT argument sequence (SUBJECT TO CHANGE) 

        # determine "idata"
        if catalog['name'].lower() == 'lasdamasgeo': 
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

    return fft_file  
