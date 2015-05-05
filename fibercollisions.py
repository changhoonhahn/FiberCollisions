import numpy as np 
import os.path
import time 
import subprocess
import cosmolopy as cosmos
import fibcol_nbar as fc_nbar
import fibcol_data as fc_data
import fibcol_spec as fc_spec
import fibcol_utility as fc_util
import fibcol_fft as fc_fft
import plot_fibcol as fc_plot

# Functions -----------------------------------------------------------------
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

def fibcoll_data_prep(DorR, silent=True, **cat_corr): 
    '''  Construct mock/random data for Pk/Bk with corrections 
    checks if data file exists, if it doesn't, makes it 

    Parameters
    ----------
    DorR : 'data' or 'random' 
    cat_corr : catalog and correction dictionary
    silent : print or not print 
    
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    if catalog['name'].lower() == 'qpm': 
        # QPM -------------------------------------------------------------------
        # QPM does not have nbar(z) correction so there's no need to append corrected nbar(z) 
        # Instead we just use the given nbar(z) 
        if DorR.lower() == 'data':      # Mock -----------------------------------
            # check that the mock file exists 
            mock_file = fc_data.get_galaxy_data_file('data', **cat_corr)

            if silent == False: 
                print mock_file 

            if os.path.isfile(mock_file) == False:  
                print 'Constructing ', mock_file
                mock = fc_data.galaxy_data('data', **cat_corr) 

        elif DorR.lower() == 'random': 
            # Random --------------------------------------------

            # corrected random file 
            corr_rand_file = fc_data.get_galaxy_data_file('random', **cat_corr) 
            
            if silent == False: 
                print corr_rand_file

            if os.path.isfile(corr_rand_file) == False: 
                print "Constructing ", corr_rand_file, ' (Will take a while!)'
                corr_rand = fc_data.galaxy_data('random', **cat_corr) 

    else:       # upweight, peak, peaktest
        if DorR.lower() == 'data':
            # Mock ---------------------------------------------

            mock_file = fc_data.get_galaxy_data_file('data', **cat_corr)    # file name 

            if silent == False: 
                print mock_file

            #if os.path.isfile(mock_file) == False:      # check if it exists 
            print 'Constructing ', mock_file
            mock = fc_data.galaxy_data('data', clobber=True, **cat_corr) 
            
            if catalog['name'].lower() == 'tilingmock': 
                # check if corrected nbar is appended
                if os.path.isfile(mock.file_name+'.corrnbar') == False:       
                    # if corrected nbar is not appended
                    # nbar does not change from peak correction so there is no need to append corrected nbar
                    print "appending corrected nbar to mock ", mock.file_name 
                    fc_nbar.append_corr_nbar('data', sanitycheck=False, **cat_corr)

        elif DorR.lower() == 'random': 
            # Random -----------------------------------------------
            # corrected random file 
            corr_rand_file = fc_data.get_galaxy_data_file('random', **cat_corr) 
            
            if silent == False: 
                print corr_rand_file

            if os.path.isfile(corr_rand_file) == False: # if there is no corrected random file

                print "Constructing ", corr_rand_file, ' (Will take a while!)'
                corr_rand = fc_data.galaxy_data('random', **cat_corr) 
            
            if catalog['name'].lower() == 'tilingmock': 
                # does append corrected nbar corrected random file exist? 
                if os.path.isfile(corr_rand.file_name+'.corrnbar') == False: 
                    # if not 
                    print "appending corrected nbar to corrected random ", corr_rand.file_name
                    fc_nbar.append_corr_nbar('random', sanitycheck=False, **cat_corr)  
# ----------------------------------------------------------------------------
# Las Damas Geo 
def lasdamasgeo_fibcoll_pk(i_mock, corr, quad=False): 
    '''
    Compute fibercollision corrected P(k)
    '''
    # set catalog_correction dictionary
    catalog = {'name':'lasdamasgeo'}
    correction = corr
    # some constants throughout the code
    spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':960, 'quad': quad} 

    cat_corr = {'catalog':catalog, 'correction': correction, 'spec': spec}
    # random data ------------------------------------------------
    fibcoll_data_prep('random', **cat_corr) 

    # build random FFT  
    rand_fft_file = fc_fft.build_fibcol_fft('random', **cat_corr)
    
    # mock data ---------------------------------------------------
    for letter in ['a', 'b', 'c', 'd']: 
        i_cat_corr = {'catalog':catalog, 'correction': correction, 'spec': spec}
        i_cat_corr['catalog']['n_mock'] = i_mock 
        i_cat_corr['catalog']['letter'] = letter

        fibcoll_data_prep('data', **i_cat_corr) 
        
        # build mock FFT
        fft_file = fc_fft.build_fibcol_fft('data', **i_cat_corr) 
        print 'Constructing ', fft_file 
        
        power_file = fc_spec.build_fibcol_power(**i_cat_corr) 
        print 'Constructing ', power_file 

        # in order to safe memory delete data FFT file which can be generated quickly
        cleanup_cmd = ''. join(['rm ', fft_file])
        print cleanup_cmd 
        os.system(cleanup_cmd) 

"""
def lasdamasgeo_fibcoll_pk_Igal_Iran(n_mocks, corr): 
    '''
    Compute fibercollision corrected P(k)
    '''
    # set catalog_correction dictionary
    catalog = {'name':'lasdamasgeo'}
    correction = {'name':corr} 
    # some constants throughout the code
    spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 
    
    if (corr.lower() == 'peak') or (corr.lower() == 'peaknbar') or (corr.lower() == 'peaktest') or (corr.lower() == 'allpeak'): 
        correction['sigma'] = 5.3
        correction['fpeak'] = 0.5           # hardcoded

    cat_corr = {'catalog':catalog, 'correction': correction, 'spec': spec}

    # mock data ---------------------------------------------------
    for i_mock in range(1, n_mocks+1): 
        for letter in ['a', 'b', 'c', 'd']: 
            i_cat_corr = {'catalog':catalog, 'correction': correction, 'spec': spec}
            i_cat_corr['catalog']['n_mock'] = i_mock 
            i_cat_corr['catalog']['letter'] = letter
            
            power_file = build_fibcol_power_Igal_Iran(**i_cat_corr) 
            print 'Constructing ', power_file 
"""

def lasdamasgeo_fibcoll_bispec(n_mocks, corr): 
    '''
    Compute fibercollision corrected P(k)
    '''
    # set catalog_correction dictionary
    catalog = {'name':'lasdamasgeo'}
    correction = {'name':corr} 
    # some constants throughout the code
    spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 
    
    if (corr.lower() == 'peak') or (corr.lower() == 'peaknbar') or (corr.lower() == 'peaktest'): 
        correction['sigma'] = 5.3
        correction['fpeak'] = 1.0           # hardcoded

    cat_corr = {'catalog':catalog, 'correction': correction, 'spec': spec}
    # random data ------------------------------------------------
    fibcoll_data_prep('random', **cat_corr) 

    # build random FFT  
    rand_fft_file = build_fibcol_fft('random', **cat_corr)

    # mock data ---------------------------------------------------
    for i_mock in range(1, n_mocks+1): 
        for letter in ['a', 'b', 'c', 'd']: 
            i_cat_corr = {'catalog':catalog, 'correction': correction, 'spec': spec}
            i_cat_corr['catalog']['n_mock'] = i_mock 
            i_cat_corr['catalog']['letter'] = letter

            fft_file = get_fibcol_fft_file('data', **i_cat_corr) 
            if os.path.isfile(fft_file): 
                # build mock FFT if it doesn't exist
                fft_file = build_fibcol_fft('data', **i_cat_corr) 
                print 'Constructing ', fft_file 
            
            bispec_file = build_fibcol_bispec(**i_cat_corr) 
            print 'Constructing ', bispec_file 

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

# ------------------------------------------------------------------------------------
# Tiling Mock
def tilingmock_fibcoll_pk(corr, quad=False): 
    '''
    Compute fibercollision corrected P(k)
    '''
    # set catalog_correction dictionary
    catalog = {'name':'tilingmock'}
    correction = corr
    # some constants throughout the code
    spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':960, 'quad':quad} 

    cat_corr = {'catalog':catalog, 'correction': correction, 'spec': spec}
    # random data ------------------------------------------------
    fibcoll_data_prep('random', silent=False, **cat_corr) 

    # build random FFT  
    rand_fft_file = fc_fft.build_fibcol_fft('random', **cat_corr)

    # mock data ---------------------------------------------------
    fibcoll_data_prep('data', **cat_corr) 
            
    # build mock FFT
    fft_file = fc_fft.build_fibcol_fft('data', **cat_corr) 
    print 'Constructing ', fft_file 
    
    power_file = fc_spec.build_fibcol_power(**cat_corr) 
    print 'Constructing ', power_file 

# QPM --------------------------------------------------------------------------------
def qpm_fibcoll_pk(i_mock, corr, quad=False): 
    '''
    Compute fibercollision corrected P(k)
    '''
    # set catalog_correction dictionary
    catalog = {'name':'qpm'}
    correction = corr
    # some constants throughout the code
    spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360, 'quad':quad}

    i_catalog = catalog.copy() 
    i_catalog['n_mock'] = i_mock 

    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec': spec}
    # random data ------------------------------------------------
    fibcoll_data_prep('random', silent=False, **i_cat_corr) 

    # build random FFT  
    rand_fft_file = fc_fft.build_fibcol_fft('random', **i_cat_corr)

    # mock data ---------------------------------------------------
    fibcoll_data_prep('data', **i_cat_corr) 
            
    # build mock FFT
    fft_file = fc_fft.build_fibcol_fft('data', **i_cat_corr) 
    print 'Constructing ', fft_file 
    
    power_file = fc_spec.build_fibcol_power(**i_cat_corr) 
    print 'Constructing ', power_file 

# PATCHY ----------------------------------------------------------------------------
def patchy_fibcoll_pk(i_mock, corr, quad=False): 
    ''' Compute fibercollision corrected P(k) for PATCHY mocks 
    '''
    # set catalog_correction dictionary
    catalog = {'name':'patchy'}
    correction = corr
    # some constants throughout the code
    spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360, 'quad':quad}

    i_catalog = catalog.copy() 
    i_catalog['n_mock'] = i_mock 

    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec': spec}
    # random data ------------------------------------------------
    fibcoll_data_prep('random', silent=False, **i_cat_corr) 

    # build random FFT  
    rand_fft_file = fc_fft.build_fibcol_fft('random', **i_cat_corr)

    # mock data ---------------------------------------------------
    fibcoll_data_prep('data', **i_cat_corr) 
            
    # build mock FFT
    fft_file = fc_fft.build_fibcol_fft('data', **i_cat_corr) 
    print 'Constructing ', fft_file 
    
    power_file = fc_spec.build_fibcol_power(**i_cat_corr) 
    print 'Constructing ', power_file 

#  Average P(k) (both monopole and quadrupole) ---------
def build_avg_Pk(n_mock, quad=False, **cat_corr):           
    ''' Calculate average P(k) using n_mocks
    '''

    catalog = cat_corr['catalog'] 
    correction = cat_corr['correction']
    
    if 'spec' not in cat_corr.keys(): 
        if catalog['name'].lower() == 'tilingmock': 
            spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360, 
                    'quad': quad} 
        else: 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360, 
                    'quad': quad} 
    else:  
        spec = cat_corr['spec']
    
    # Compute total( P(k) )
    if catalog['name'].lower() == 'qpm': 
        # QPM 
        for i_mock in range(1, n_mock+1): 
            i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
            (i_cat_corr['catalog'])['n_mock'] = i_mock
                
            power_i = fc_spec.Spec('power', **i_cat_corr)
            power_i.readfile()
            
            if spec['quad'] == False: 
                Pk_i = power_i.Pk
            else: 
                Pk_i = power_i.P2k

            try: 
                tot_Pk
            except NameError: 
                tot_Pk = Pk_i
            else: 
                tot_Pk = tot_Pk + Pk_i

        n_file = n_mock 

    elif catalog['name'].lower() == 'lasdamasgeo': 
        # LasDamasGeo

        n_file = 0
        for i_mock in range(1, n_mock+1): 
            for letter in ['a', 'b', 'c', 'd']: 
                i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
                (i_cat_corr['catalog'])['n_mock'] = i_mock
                (i_cat_corr['catalog'])['letter'] = letter 
                
                power_i = fc_spec.Spec('power', **i_cat_corr)
                power_i.readfile()
            
                if spec['quad'] == False: 
                    Pk_i = power_i.Pk
                else: 
                    Pk_i = power_i.P2k

                try: 
                    tot_Pk
                except NameError: 
                    tot_Pk = Pk_i
                else: 
                    tot_Pk = tot_Pk + Pk_i
                n_file += 1 

    elif catalog['name'].lower() == 'tilingmock': 
        n_file = 1
        i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec': spec}
        power_i = fc_spec.Spec('power', **i_cat_corr)
        power_i.readfile()
            
        if spec['quad'] == False: 
            tot_Pk = power_i.Pk
        else: 
            tot_Pk = power_i.P2k

    else: 
        raise NotImplementedError('errorerereor') 

    avg_k = power_i.k  
    avg_Pk = tot_Pk/np.float(n_file)
    
    # output to file 
    avg_file_name = avg_Pk_file(n_mock, **cat_corr)     
    np.savetxt(avg_file_name, np.c_[avg_k, avg_Pk], 
            fmt=['%10.5e', '%10.5e'], delimiter='\t') 

def avg_Pk_file(n_mock, quad=False, **cat_corr):            # avg P(k) file name 
    ''' Return avareage P(k) file name (both monopole and quadrupole) 

    Paramters
    ---------
    n_mock : number of mocks
    cat_corr : catalog correction dictionary 

    '''
    catalog = cat_corr['catalog'] 

    if 'spec' not in cat_corr.keys(): 
        if catalog['name'].lower() == 'qpm': 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360, 
                    'quad': quad} 
            (cat_corr['catalog'])['n_mock'] = 1
        elif catalog['name'].lower() == 'lasdamasgeo': 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360, 
                    'quad': quad} 
            (cat_corr['catalog'])['n_mock'] = 1 
            (cat_corr['catalog'])['letter'] = 'a'
        elif catalog['name'].lower() == 'tilingmock': 
            spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360, 
                    'quad': quad} 
        else:
            raise NameError("error error") 
        cat_corr['spec'] = spec
    else: 
        if catalog['name'].lower() == 'qpm': 
            (cat_corr['catalog'])['n_mock'] = 1

        elif catalog['name'].lower() == 'lasdamasgeo': 
            (cat_corr['catalog'])['n_mock'] = 1 
            (cat_corr['catalog'])['letter'] = 'a'

    power_i = fc_spec.Spec('power', **cat_corr)
    
    #a0.6452_0045.dr12d_cmass_ngc.vetoed.fibcoll.gauss.peaknbar.sigma4.38.fpeak1.0.dat
    indiv_filename = power_i.file_name           

    if catalog['name'].lower() == 'qpm': 
        #a0.6452_0045.dr12d_cmass_ngc.vetoed.fibcoll.gauss.peaknbar.sigma4.38.fpeak1.0.dat
        indiv_filename = power_i.file_name           

        # get average P(k) file name 
        avg_file_name = indiv_filename.split('6452_')[0]+'6452_'+str(n_mock)+'mockavg.dr12d'+(indiv_filename.split('6452_')[1]).split('.dr12d')[1]

    elif catalog['name'].lower() == 'lasdamasgeo': 

        indiv_filename = power_i.file_name          # sdssmock_gamma_lrgFull_zm_oriana19b_no.rdcz.fibcoll.dat.expon.peakshot.sigma7.0.fpeak0.8.corrnbar 

        avg_file_name = indiv_filename.split('oriana')[0]+'oriana'+str(n_mock)+'mockavg_no.'+(indiv_filename.split('oriana')[1]).split('_no.')[1]

    elif catalog['name'].lower() == 'tilingmock': 

        avg_file_name = power_i.file_name          # just the file itself since it only has 1 file

    return avg_file_name

def avgPk(n_mock, clobber=False, **cat_corr):
    ''' Wrapper that returns average(k) and average(P(k))
    '''

    avg_file_name = avg_Pk_file(n_mock, **cat_corr) 

    i_cat_corr = cat_corr.copy()                # get the highest number mock power 
    (i_cat_corr['catalog'])['n_mock'] = n_mock  
    power_hii = fc_spec.Spec('power', **i_cat_corr)

    # check that the file exists and check that the file 
    # is more updated that the highest number mock power file 
    high_file_mod_time = os.path.getmtime(power_hii.file_name)
    if os.path.isfile(avg_file_name) == False: 
        avg_file_mod_time = 0 
    else: 
        avg_file_mod_time = os.path.getmtime(avg_file_name)

    # if code was changed since exe file was last compiled then compile fft code 
    if (avg_file_mod_time < high_file_mod_time) or (clobber == True): 
        build_avg_Pk(n_mock, **cat_corr) 
    
    avg_k, avg_Pk = np.loadtxt(avg_file_name, unpack=True, usecols=[0, 1]) 

    return [avg_k, avg_Pk]

def avgP_interp(k, avg_k, avg_Pk):          # get average P(k) of n_mocks at k using numpy interpolate  (ujust for convenience)
    k_avg_Pk = np.interp(k, avg_k, avg_Pk)
    return k_avg_Pk

#  DELTA PK ----------
def deltaP_file(n_mock, **cat_corr):        
    ''' Return name of delta P(k) fil 
    '''
    catalog = cat_corr['catalog']
    spec = cat_corr['spec']

    avg_power_file_name = avg_Pk_file(n_mock, **cat_corr) 
    
    if spec['quad'] == True: 
        delP_str = 'delP2_'
        quad_str = 'quad_'
    else: 
        delP_str = 'delP_'
        quad_str = ''


    if catalog['name'].lower() == 'qpm':    # QPM 
        Pk_var_file_name = (delP_str+'a0').join(avg_power_file_name.split('power_'+quad_str+'a0')) 
    elif catalog['name'].lower() == 'lasdamasgeo':  # LasDamasGeo
        Pk_var_file_name = (delP_str+'sdssmock').join(avg_power_file_name.split('power_'+quad_str+'sdssmock')) 
    elif catalog['name'].lower() == 'tilingmock': 
        Pk_var_file_name = (delP_str+'cmass').join(avg_power_file_name.split('power_'+quad_str+'cmass')) 

    return Pk_var_file_name

def build_deltaP(n_mock, quad=False, **cat_corr):             
    '''Calculated Delta P for a given correction method
    '''
    catalog = cat_corr['catalog']
    
    # get average P(k) file name 
    avg_power = avgPk(n_mock, **cat_corr) 
    avg_k = avg_power[0]
    avg_Pk = avg_power[1]

    if 'spec' not in cat_corr.keys(): 
        if catalog['name'].lower() == 'qpm': 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360, 
                    'quad': quad} 
        elif catalog['name'].lower() == 'lasdamasgeo': 
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360, 
                    'quad': quad}  
        elif catalog['name'].lower() == 'tilingmock':
            spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360, 
                    'quad': quad} 
        else: 
            raise NotImplementedError('asdfasdf') 
    else: 
        spec = cat_corr['spec']

    # Calculate DeltaP for QPM 
    if catalog['name'].lower() == 'qpm': 

        for i_mock in range(1, n_mock+1): 
            i_cat_corr = cat_corr.copy() 
            (i_cat_corr['catalog'])['n_mock'] = i_mock 
            i_cat_corr['spec'] = spec

            power_i = fc_spec.Spec('power', **i_cat_corr)
            power_i.readfile()
            
            if spec['quad'] == False: 
                Pk_i = power_i.Pk
            else: 
                Pk_i = power_i.P2k

            try: 
                Pk_var
            except NameError: 
                Pk_var = (avg_Pk - Pk_i)**2
            else: 
                Pk_var = Pk_var + (avg_Pk - Pk_i)**2

        n_file = n_mock

    elif catalog['name'].lower() == 'lasdamasgeo': 

        n_file = 0 
        for i_mock in range(1, n_mock+1): 
            for letter in ['a', 'b', 'c', 'd']: 
                # set individual i_Cat_corr
                i_cat_corr = cat_corr.copy() 
                (i_cat_corr['catalog'])['n_mock'] = i_mock 
                (i_cat_corr['catalog'])['letter'] = letter 
                i_cat_corr['spec'] = spec

                power_i = fc_spec.Spec('power', **i_cat_corr)
                power_i.readfile()
                
                if spec['quad'] == False: 
                    Pk_i = power_i.Pk
                else: 
                    Pk_i = power_i.P2k

                try: 
                    Pk_var
                except NameError: 
                    Pk_var = (avg_Pk - Pk_i)**2
                else: 
                    Pk_var = Pk_var + (avg_Pk - Pk_i)**2
                n_file = n_file+1
    elif catalog['name'].lower() == 'tilingmock': 
        # TM deltaP is QPMdeltaP scaled by the ratio of the true power spectrum 
        qpm_n_mock = 100
        qpm_cat_corr = {'catalog': {'name': 'qpm'}, 
                'correction': {'name': 'true'}, 
                'spec': spec} 
        qpm_Pk = avgPk(qpm_n_mock, **qpm_cat_corr)
        deltaPk = deltaP(qpm_n_mock, **qpm_cat_corr) 

        if (cat_corr['correction'])['name'].lower() == 'true': 
            tm_true_Pk = avg_power
        else: 
            tm_true_cat_corr = cat_corr.copy()  
            tm_true_cat_corr['correction'] = {'name': 'true'}
            tm_true_Pk = avgPk(1, **tm_true_cat_corr) 
        
        # QPM deltaP for TM avg_k 
        qpm_delP = np.array([deltaP_interp((tm_true_Pk[0])[i], deltaPk[0], deltaPk[1]) for i in range(len(tm_true_Pk[0]))])

        Pk_var = qpm_delP*(tm_true_Pk[1]/qpm_Pk[1])     # square root already taken into account from QPM DeltaP!
        
        n_file = 1
    else: 
        raise NameError("only QPM and LasDamasGeo have variance") 
    
    if catalog['name'].lower() != 'tilingmock':         # if not tiling mock then square root it
        Pk_var = np.sqrt(Pk_var/np.float(n_file))
    
    # output  average P(k) 
    Pk_var_file_name = deltaP_file(n_mock, **cat_corr)
    np.savetxt(Pk_var_file_name, np.c_[avg_k, Pk_var], fmt=['%10.5e', '%10.5e'], delimiter='\t') 

def deltaP(n_mock, clobber=False, **cat_corr): 
    ''' Calculate delta P for n_mocks for given catalog and correction scheme 

    Parameters
    ----------
    n_mock : # of mocks 
    clobber : If True then recalculates deltaP 
    cat_corr : catalog correction dictionary 

    '''

    avgP_file_name = avg_Pk_file(n_mock, **cat_corr)
    deltaP_file_name = deltaP_file(n_mock, **cat_corr)  
    
    if os.path.isfile(avgP_file_name) == True: 
        avgP_file_mod_time = os.path.getmtime(avgP_file_name)
    else: 
        build_avg_Pk(n_mock, **cat_corr)
        avgP_file_mod_time = os.path.getmtime(avgP_file_name)

    # check that the file exists and check that the file is more updated that the highest number mock power file 
    if (os.path.isfile(deltaP_file_name) == False) or \
            (clobber == True): 
        deltaP_file_mod_time = 0 

    else: 
        deltaP_file_mod_time = os.path.getmtime(deltaP_file_name)
    
    if (deltaP_file_mod_time < avgP_file_mod_time) or (clobber == True): 
        build_deltaP(n_mock, **cat_corr) 

    avg_k, delta_P = np.loadtxt(deltaP_file_name, unpack=True, usecols=[0,1])
    
    return [avg_k, delta_P]

def deltaP_interp(k, avg_k, Pk_var):       # given avg_k and DeltaP, inteerpolate to k 
    k_Pk_var = np.interp(k, avg_k, Pk_var)

    return k_Pk_var

#-----------------------------------------------------------------------------
def pk_fibcol_comp(catalog_name, n_mock, corr_methods): 
    '''
    Comparison of P(k)avg for different fibercollision correciton methods. Calculate minimum chisquared correction method 
    Don't include True
    '''
    catalog = {'name': catalog_name}    # catalog dict
    # hardcoded default power/bispectrum box settings ----------------------------------------------------------------
    if catalog_name.lower() in ('lasdamasgeo', 'qpm'): 
        spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 
    elif catalog_name.lower() == 'tilingmock': 
        spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360} 
    else: 
        raise NameError('not coded yet') 
    #----------------------------------------------------------------------------------------------------------------
   
    # Import average True P(K) 
    true_cat_corr = {}
    true_cat_corr['catalog'] = catalog
    true_cat_corr['correction'] = {'name': 'true'}
    true_cat_corr['spec'] = spec
    avg_Pk_true_file = avg_Pk_file(n_mock, **true_cat_corr) 
    avg_k_true, avg_Pk_true = np.loadtxt(avg_Pk_true_file, unpack=True, usecols = [0,1]) 

    # Import P(k) variance for true
    if catalog['name'].lower() == 'qpm': 
        Pk_var_file_name = 'delP_a0'.join(avg_Pk_true_file.split('power_a0')) 
    elif catalog['name'].lower() == 'lasdamasgeo': 
        Pk_var_file_name = 'delP_sdssmock'.join(avg_file_name.split('power_sdssmock')) 
    elif catalog['name'].lower() == 'tilingmock':  
        pass
    else: 
        raise NameError("asdfasdfasdfasdfas") 

    if catalog['name'].lower != 'tilingmock': 
        Pk_var = np.loadtxt(Pk_var_file_name, unpack=True, usecols=[1])
    else: 
        qpm_cat_corr = {'catalog': {'name': 'qpm'}, 'correction': {'name':'true'}, 'spec': spec} 
        qpm_Pk = np.array([
            avgP_interp(avg_k_true[i], 45, **qpm_cat_corr) for i in range(len(avg_k_true))
            ])
        qpm_delP = np.array([
            deltaP_interp(avg_k_true[i], 45, **qpm_cat_corr) for i in range(len(avg_k_true))
            ])
        Pk_var = (qpm_delP/qpm_Pk) * avg_Pk_true

    chi2 = np.zeros(len(corr_methods)-1)
    for i_corr, correction in enumerate(corr_methods):     # loop through correction methods
        # for LasDamasGeo 
        if catalog_name.lower() == 'lasdamasgeo': 
            n_file = 0
            for i_mock in range(1,n_mock+1):                       # compute average[P(k)] for each correction method
                for letter in ['a', 'b', 'c', 'd']: 
                    # set catalog correction dictionary for specific file 
                    i_catalog = catalog.copy() 
                    i_catalog['n_mock'] = i_mock
                    i_catalog['letter'] = letter
                    i_cat_corr = {'catalog': i_catalog, 'correction': correction, 'spec':spec}
                    
                    if correction['name'] != 'Upweight Igal_Irand': 
                        # import power spectrum 
                        power_i = fc_spec.Spec('power', **i_cat_corr)
                    else: 
                        i_cat_corr = {'catalog': i_catalog, 'correction': {'name':'upweight'}, 'spec':spec}
                        power_i = fc_spec.Spec('power', Igal_Irand==True, **i_cat_corr)
                    power_i.readfile()
                    
                    try: 
                        avg_k 
                    except NameError: 
                        avg_k = power_i.k
                        sum_Pk = power_i.Pk
                    else: 
                        sum_Pk = sum_Pk + power_i.Pk

                    n_file = n_file+1

        # Tiling Mocks
        elif catalog_name.lower() == 'tilingmock': 
            i_cat_corr = {'catalog': catalog, 'correction': correction, 'spec':spec}
            # read tiling mock P(k)
            power = fc_spec.Spec('power', **i_cat_corr)
            power.readfile()

            avg_k = power.k
            sum_Pk = power.Pk
            n_file = 1          # only one file 

        avg_Pk = np.array([sum_Pk[i]/np.float(n_file) for i in range(len(sum_Pk))])    # average P(k)
        
        chi2[i_corr-1] = np.sum((avg_Pk - avg_Pk_true)**2/(Pk_var)**2) 

        del avg_k
        del avg_Pk

    min_index = np.argmin(chi2)

    print corr_methods[min_index+1]

'''
def pk_ldg_peakshot_expon_grid(): 
    spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360}
    
    ldg_sigmas = [5.0+0.1*np.float(i) for i in range(0,21)]
    ldg_fpeaks = [0.3+0.1*np.float(i) for i in range(0,6)]

    for ldg_sigma in ldg_sigmas: 
        for ldg_fpeak in ldg_fpeaks: 
            print 'SIGMA = ', ldg_sigma, ' FPEAK = ', ldg_fpeak
            corr = {'name': 'peakshot', 'sigma':ldg_sigma, 'fpeak':ldg_fpeak, 'fit':'expon'}
            #lasdamasgeo_fibcoll_pk(40, corr) 
            for i_mock in range(1,41): 
                for letter in ['a', 'b', 'c', 'd']: 
                    cat = {'name': 'lasdamasgeo', 'n_mock': i_mock, 'letter': letter} 
                    cat_corr = {'catalog': cat, 'correction': corr, 'spec': spec} 
                    power = fc_spec.Spec('power', **cat_corr) 
                    power_file = power.file_name 

                    if os.path.isfile(power_file) == False: 
                        print 'Constructign ', power_file 
                        lasdamasgeo_fibcoll_pk(i_mock, corr)

            # compute the average for the correction: 
            cat = {'name': 'lasdamasgeo'}
            cat_corr = {'catalog': cat, 'correction': corr, 'spec': spec}  
            build_avg_Pk(40, **cat_corr)
   
    # import average true P(k)
    true_cat_corr = {'catalog': {'name': 'lasdamasgeo'}, 'correction': {'name': 'true'}, 'spec': spec} 
    avg_Pk_true_file = avg_Pk_file(40, **true_cat_corr) 
    if os.path.isfile(avg_Pk_true_file) == False:
        print 'Constructing ', avg_Pk_true_file
        build_avg_Pk(40, **true_cat_corr) 

    avg_Pk_true = np.loadtxt(avg_Pk_true_file, unpack=True, usecols=[1])
    
    # import delta P(k) for true
    Pk_var_file_name = 'delP_sdssmock'.join(avg_Pk_true_file.split('power_sdssmock')) 
    if os.path.isfile(Pk_var_file_name) == False: 
        deltaP(40, **true_cat_corr)

    delPk = np.loadtxt(Pk_var_file_name, unpack=True, usecols=[1]) 

    chi2 = [] 
    chi2_param = [] 
    for ldg_sigma in ldg_sigmas: 
        for ldg_fpeak in ldg_fpeaks: 
            chi2_param.append([ldg_sigma, ldg_fpeak])
            print 'SIGMA = ', ldg_sigma, ' FPEAK = ', ldg_fpeak
            cat = {'name': 'lasdamasgeo'}
            corr = {'name': 'peakshot', 'sigma':ldg_sigma, 'fpeak':ldg_fpeak, 'fit':'expon'}
            cat_corr = {'catalog': cat, 'correction': corr, 'spec': spec}  

            avg_Pk_filename = avg_Pk_file(40, **cat_corr)

            avg_k, avg_Pk = np.loadtxt(avg_Pk_filename, unpack=True, usecols=[0,1])

            chi2.append(np.sum((avg_Pk_true-avg_Pk)**2/(delPk**2)))
    
    min_chi2_index = chi2.index(min(chi2))
    print 'best fit parameters ------- '
    print 'sigma = ', chi2_param[min_chi2_index][0]
    print 'fpeak = ', chi2_param[min_chi2_index][1]

def pk_tm_peakshot_expon_grid(): 
    qpm_nmock = 45  # hardcoded value 

    spec = {'P0': 20000, 'sscale':4000.0, 'Rbox':2000.0, 'box':4000, 'grid':360} 

    tm_sigmas = [5.0+0.1*np.float(i) for i in range(0,21)]
    tm_fpeaks = [0.3+0.1*np.float(i) for i in range(0,6)]
    for tm_sigma in tm_sigmas: 
        for tm_fpeak in tm_fpeaks: 
            print 'SIGMA = ', tm_sigma, ' FPEAK = ', tm_fpeak
            corr = {'name': 'peakshot', 'sigma':tm_sigma, 'fpeak':tm_fpeak, 'fit':'expon'}
            cat_corr = {'catalog': {'name':'tilingmock'}, 'correction': corr, 'spec': spec} 

            power = fc_spec.Spec('power', **cat_corr) 
            power_file = power.file_name 

            if os.path.isfile(power_file) == False: 
                print 'Constructign ', power_file 
                tilingmock_fibcoll_pk(corr) 
    
    # calculate which one has the best fit 
    cat_corr = {'catalog': {'name':'tilingmock'}, 'correction': {'name': 'true'}, 'spec': spec} 
    power_true = fc_spec.Spec('power', **cat_corr) 
    power_true.readfile()
    k = power_true.k
    pk_true = power_true.Pk

    # Scaled DeltaP/P using QPM mocks
    qpm_spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':360} 
    
    # First get average QPM Pk_true
    qpm_tot_Pk_true = np.zeros(len(k)) 
    for i_mock in range(1, qpm_nmock+1): 
        i_cat_corr = {'catalog': {'name':'qpm', 'n_mock': i_mock}, 'correction': {'name':'true'}, 'spec':qpm_spec}
            
        power_i = fc_spec.Spec('power', **i_cat_corr)
        power_i.readfile()
        
        Pk_i = power_i.Pk
        qpm_tot_Pk_true = qpm_tot_Pk_true + Pk_i

    qpm_avg_Pk_true = qpm_tot_Pk_true/np.float(qpm_nmock)
    
    # Calculate DeltaP for QPM 
    qpm_Pk_var = np.zeros(len(k)) 
    for i_mock in range(1, qpm_nmock+1): 
        i_cat_corr = {'catalog': {'name':'qpm', 'n_mock': i_mock}, 'correction': {'name':'true'}, 'spec':qpm_spec}
        power_i = fc_spec.Spec('power', **i_cat_corr)
        power_i.readfile()
        
        Pk_i = power_i.Pk
        qpm_Pk_var = qpm_Pk_var + (qpm_avg_Pk_true - Pk_i)**2

    qpm_Pk_var = np.sqrt(qpm_Pk_var/np.float(qpm_nmock))
    tm_Pk_var = (qpm_Pk_var/qpm_avg_Pk_true) * pk_true

    chi2 = [] 
    chi2_param = [] 
    for tm_sigma in tm_sigmas:
        for tm_fpeak in tm_fpeaks: 
            chi2_param.append([tm_sigma, tm_fpeak])
            corr = {'name': 'peakshot', 'sigma':tm_sigma, 'fpeak':tm_fpeak, 'fit':'expon'}
            cat_corr = {'catalog': {'name':'tilingmock'}, 'correction': corr, 'spec': spec} 

            power = fc_spec.Spec('power', **cat_corr) 
            power.readfile()

            k = power.k 
            pk = power.Pk

            print np.sum((pk_true/pk-1.0)**2)
            chi2.append(np.sum((pk_true-pk)**2/(tm_Pk_var**2)))
    
    min_chi2_index = chi2.index(min(chi2))
    print 'best fit parameters ------- '
    print 'sigma = ', chi2_param[min_chi2_index][0]
    print 'fpeak = ', chi2_param[min_chi2_index][1]
    corr = {'name': 'peakshot', 'sigma':chi2_param[min_chi2_index][0],'fpeak':chi2_param[min_chi2_index][1], 'fit':'expon'}
    cat_corr = {'catalog': {'name':'tilingmock'}, 'correction': corr, 'spec': spec} 

    power = fc_spec.Spec('power', **cat_corr) 

    k, pk = np.loadtxt(power.file_name, unpack=True, usecols=[0,1])

    print pk_true/pk-1.0
'''

def build_qpm():
    ''' Function to build QPM stuff
    ''' 
    pass 
    '''
    #corr = {'name': 'peaknbar', 'sigma':4.38, 'fpeak':1.0, 'fit':'gauss'}
    #corr = {'name': 'vlospeakshot', 'sigma':650, 'fpeak':0.7, 'fit':'expon'}
    #[{'name': 'true'}, {'name':'upweight'}, {'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'gauss'}]#[{'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'gauss'}]
    #[{'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'gauss'}]#, {'name': 'peakshot', 'sigma': 4.8, 'fpeak': 0.65, 'fit': 'expon'}]
            #{'name': 'true'}, {'name':'upweight'}, {'name': 'peakshot', 'sigma': 4.3, 'fpeak': 0.6, 'fit': 'gauss'}]
            #{'name': 'true'}, {'name':'upweight'}, {'name': 'floriansn'}, {'name': 'hectorsn'}] 
            #{'name': 'peakshot', 'sigma':5.0, 'fpeak':0.65, 'fit':'expon'}, {'name': 'vlospeakshot', 'sigma':580, 'fpeak':0.62, 'fit':'gauss'}] 
    #corrections = [{'name': 'true'}, {'name': 'upweight'}, {'name': 'vlospeakshot', 'sigma':580, 'fpeak':0.62, 'fit':'gauss'}, {'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'gauss'}]
    #corrections = [{'name': 'peakshot', 'fit': 'true', 'fpeak':0.65}]
    ##corrections = [{'name': 'peakshot', 'fit': 'true', 'fpeak':1.0}]
    #corrections = [{'name': 'tailupw'}] 
    #####corrections = [ 
    #####        {'name': 'true'}, 
    #####        {'name': 'upweight'}, 
    #####        {'name': 'peakshot', 'fit': 'gauss', 'sigma': 4.4, 'fpeak':0.65}, 
    #####        {'name': 'floriansn'}, 
    #####        {'name': 'hectorsn'}
    #####        ]
    #####for i_mock in range(51, 101): 
    #####    for corr in corrections: 
    #####        qpm_fibcoll_pk(i_mock, corr, quad=False) 
   
    #####corrections = [{'name': 'peaknbar', 'sigma': 4.4, 'fpeak': 1.0, 'fit': 'gauss'}, {'name': 'shotnoise'}] 
    #####
    #####for i_mock in range(1, 101): 
    #####    for corr in corrections: 
    #####        qpm_fibcoll_pk(i_mock, corr, quad=False) 
    '''

def build_patchy(): 
    corrs = [{'name': 'true'}, {'name': 'upweight'}] 
    for i_mock in np.arange(1,11): 
        for corr in corrs: 
            patchy_fibcoll_pk(i_mock, corr, quad=False) 

if __name__=='__main__': 
    build_patchy()

    #for i_mock in np.arange(1,2): 
    #    fc_data.galaxy_data('data', readdata=True, clobber=True, **cat_corr) 

    #pk_ldg_peakshot_expon_grid()
    #qpm_avgP(45, {'name':'true'})
    #qpm_deltaP(45, {'name':'true'})
    #pk_tm_peakshot_expon_grid()
    
    #corr = {'name': 'upweight'}
    #qpm_fibcoll_pk(2, corr) 
    #cat_corr = {'catalog': {'name':'lasdamasgeo'}, 'correction': corr} 
    #avg_Pk(40, **cat_corr)
    #corr = {'name': 'peaknbar', 'sigma':6.9, 'fpeak':1.0, 'fit':'expon'}
    #lasdamasgeo_fibcoll_pk(40, corr) 

    ############################################################################################# QPM 

    ############################################################################################# LasDamasGeo  
    #cat_corr = {'catalog': {'name': 'qpm'}, 'correction': {'name': 'true'}} 
    #build_avg_Pk(100, **cat_corr)
    #deltaP(100, **cat_corr)
    #[{'name': 'vlospeakshot', 'fpeak': 0.6, 'fit':'real'}, {'name': 'vlospeakshot', 'fpeak': 0.7, 'fit':'real'}, {'name': 'vlospeakshot', 'fpeak': 0.8, 'fit':'real'}] #{'name': 'hectorsn'}, {'name': 'floriansn'}, 
    #corrections = [{'name': 'peaknbar', 'sigma': 6.5, 'fpeak': 1.0, 'fit': 'gauss'}, {'name': 'shotnoise'}]
    #{'name': 'floriansn'}, {'name': 'hectorsn'}, {'name': 'true'}, {'name': 'upweight'}, 

    #corrections = [{'name': 'true'}, {'name': 'upweight'}, {'name': 'peakshot', 'sigma': 6.5, 'fpeak': 0.76, 'fit': 'gauss'}, {'name': 'floriansn'}, {'name': 'hectorsn'}]
    #for i_mock in range(1,41): 
    #    for corr in corrections: 
    #        lasdamasgeo_fibcoll_pk(i_mock, corr, quad=False) 
    #corr = {'name': 'peakshot', 'sigma':6.9, 'fpeak':0.7, 'fit':'expon'}
    #lasdamasgeo_fibcoll_pk(40, corr) 
    #corr = {'name': 'peakshot', 'sigma':5.9, 'fpeak':0.7, 'fit':'gauss'}
    #lasdamasgeo_fibcoll_pk(1, corr) 
    
    ############################################################################################# Tiling Mock 
    #corr = {'name': 'peakshot', 'sigma':4.6, 'fpeak':0.53, 'fit':'gauss'}
    #tilingmock_fibcoll_pk(corr) 
    #corr = {'name': 'vlospeakshot', 'sigma':700, 'fpeak':0.65, 'fit':'expon'}
    #tilingmock_fibcoll_pk(corr) 
    #corrections = [{'name': 'shotnoise'}, {'name': 'peaknbar', 'sigma':4.8, 'fpeak':1.0, 'fit':'gauss'}]   #[{'name': 'peakshot', 'sigma':4.8, 'fpeak':0.63, 'fit':'gauss'}]#{'name': 'true'}
    #[{'name': 'floriansn'}, {'name': 'hectorsn'}]
    #{'name': 'upweight'}, {'name': 'peakshot', 'sigma':4.8, 'fpeak':0.62, 'fit':'gauss'}]#{'name': 'true'}
    #corrections = [{'name': 'peakshot', 'sigma':4.8, 'fpeak':0.63, 'fit':'gauss'}]
    #corrections = [{'name': 'true'}, {'name': 'upweight'}, {'name': 'peakshot', 'sigma':4.8, 'fpeak':0.0, 'fit':'gauss'}, {'name': 'floriansn'}, {'name': 'hectorsn'}]
    ##for corr in corrections: 
    #    tilingmock_fibcoll_pk(corr, quad=False)
