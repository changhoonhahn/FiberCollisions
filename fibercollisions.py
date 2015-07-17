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

def fibcoll_data_prep(DorR, silent=True, clobber=False, **cat_corr): 
    '''  Construct mock/random data for Pk/Bk with corrections.
    Code checks if data file exists, if it doesn't, makes it.

    Parameters
    ----------
    DorR : 'data' or 'random' 
    cat_corr : catalog and correction dictionary
    silent : print or not print 
    
    '''
    catalog = cat_corr['catalog']
    correction = cat_corr['correction']
    
    if DorR.lower() == 'data':              # Mock -----------------------------------

        mock_file = fc_data.get_galaxy_data_file('data', **cat_corr)

        if silent == False: 
            print mock_file 

        if (os.path.isfile(mock_file) == False) or (clobber == True):  
            print 'Constructing ', mock_file
            mock = fc_data.galaxy_data('data', clobber=True, **cat_corr) 
            
        if catalog['name'].lower() == 'tilingmock':     # for tiling mock 

            # check if corrected nbar is appended
            if (os.path.isfile(mock.file_name+'.corrnbar') == False) or (clobber == True):
                # if corrected nbar is not appended
                # nbar does not change from peak correction so there is no need to append corrected nbar
                print "appending corrected nbar to mock ", mock.file_name 
                fc_nbar.append_corr_nbar('data', sanitycheck=False, **cat_corr)

    elif DorR.lower() == 'random':          # Random --------------------------------------------

        # corrected random file 
        corr_rand_file = fc_data.get_galaxy_data_file('random', **cat_corr) 
        
        if not silent: 
            print corr_rand_file

        if not os.path.isfile(corr_rand_file) or clobber: 
            print "Constructing ", corr_rand_file, ' (Will take a while!)'
            corr_rand = fc_data.galaxy_data('random', clobber=True, **cat_corr) 
        
        if (catalog['name'].lower() == 'tilingmock') and (clobber == True): 
            # does append corrected nbar corrected random file exist? 
            if os.path.isfile(corr_rand.file_name+'.corrnbar') == False: 
                # if not 
                print "appending corrected nbar to corrected random ", corr_rand.file_name
                fc_nbar.append_corr_nbar('random', sanitycheck=False, **cat_corr)  
    else: 
        raise NameError("Only 'data' or 'random'")

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

def build_fibcol_pk(cat, i_mock, corr, quad=False, clobber=False, **kwargs): 
    ''' Construct P(k) (monopole or quadrupole) for a specified mock catalog 

    Parameters
    ----------
    cat : Catalog ('qpm', 'lasdamasgeo', 'tilingmock', 'nseries') 
    i_mock : Realization # of the mocks
    corr : Fiber collision correction method with specificiation 
    quad : monopole or quadrupole 

    '''
    mock_list = [] 
    if cat['name'].lower() in ('lasdamasgeo', 'ldgdownnz'):     # lasdamasgeo mocks 
        for letter in ['a', 'b', 'c', 'd']: 
            mock_list.append({'n_mock': i_mock, 'letter': letter}) 
    else: 
        mock_list.append({'n_mock': i_mock}) 
    
    # cat corr dictionary
    correction = corr   
    catalog = cat 

    if 'spec' in kwargs.keys(): 
        spec = kwargs['spec']
    else: 
        spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 
                'grid':360, 'quad':quad}

    cat_corr = {'catalog': catalog, 'correction': correction, 'spec': spec}

    # build random data 
    fibcoll_data_prep('random', silent=False, clobber=False, **cat_corr) 
    # build random FFT  
    rand_fft_file = fc_fft.build_fibcol_fft('random', **cat_corr)

    for mock_dict in mock_list: # loop through 
        # set up catalog correction dictionary 
        (cat_corr['catalog'])['n_mock'] = mock_dict['n_mock'] 
        try: 
            (cat_corr['catalog'])['letter'] = mock_dict['letter']
        except (NameError, KeyError): 
            pass 

        # mock data ---------------------------------------------------
        fibcoll_data_prep('data', clobber=clobber, **cat_corr) 
                
        # build mock FFT
        fft_file = fc_fft.build_fibcol_fft('data', **cat_corr) 
        print 'Constructing ', fft_file 
        
        power_file = fc_spec.build_fibcol_power(**cat_corr) 
        print 'Constructing ', power_file 
        
        # in order to safe memory delete data FFT file which can be generated quickly
        cleanup_cmd = ''. join(['rm ', fft_file])
        #print cleanup_cmd 
        #os.system(cleanup_cmd) 
    return 

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
    if catalog['name'].lower() in ('qpm', 'nseries'): 
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
        if catalog['name'].lower() == 'nseries': 
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
        if catalog['name'].lower() in ('qpm', 'nseries'): 
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

    elif catalog['name'].lower() == 'nseries':

        # CutskyN1.fibcoll.gauss.peakshot.sigma4.0.fpeak0.7_fidcosmo.dat
        indiv_filename = power_i.file_name           

        # get average P(k) file name 
        avg_file_name = indiv_filename.split('CutskyN')[0]+'CutskyN'+str(n_mock)+'mockavg.'+'.'.join(indiv_filename.split('.')[1:])

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
    elif catalog['name'].lower() == 'nseries':  # N series
        Pk_var_file_name = (delP_str+'Cut').join(avg_power_file_name.split('power_'+quad_str+'Cut')) 
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
    if catalog['name'].lower() in ('qpm', 'nseries'): 

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
    
    if os.path.isfile(avgP_file_name): 
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

def build_pk(catalog, n_mocks, quad=False, clobber=True, **kwargs): 
    ''' 
    Wrapper to build powerspectrum (monopole or quadrupole) for mock catalogs (QPM, Nseries, LasDamasGeo, TilingMock) 

    Parameters
    ----------
    catalog : name of catalog  
    n_mocks : # of mocks 

    Notes
    -----
    * QPM peakshot bestfit parameters: {'name': 'peakshot', 'sigma': 4.4, 'fpeak': 0.65, 'fit': 'gauss'}
    * LasDamasGeo peakshot bestfit parameters: {'name': 'peakshot', 'sigma': 6.5, 'fpeak': 0.76, 'fit': 'gauss'}
    * Nseries peakshot bestfit parameters: {'name': 'peakshot', 'sigma', 4.0, 'fpeak': 0.7, 'fit': 'gauss'}

    '''
    try: 
        Ngrid = kwargs['grid']
    except KeyError: 
        Ngrid = 360

    try:        # fiducial cosmology unless specified
        cosmology = kwargs['cosmology']
    except KeyError: 
        cosmology = 'fiducial'

    cat = {'name': catalog, 'cosmology': cosmology} 
    #corrections = [{'name': 'bigfc_peakshot', 'sigma': 6.5, 'fpeak': 0.76, 'fit': 'gauss'}]
    #corrections = [{'name': 'bigfc'}]
    corrections = [{'name': 'true'}]
    if catalog == 'cmass': 
        corrections = [{'name': 'upweight'}]

    spec = {'P0': 20000, 'sscale':3600.0, 'Rbox':1800.0, 'box':3600, 'grid':Ngrid, 'quad':quad}

    for i_mock in range(1, n_mocks+1): 
        for corr in corrections: 
            build_fibcol_pk(cat, i_mock, corr, spec=spec, clobber=clobber) 

# Compare fiber collisions of data 
def fibcol_now_fraction(**cat_corr): 
    ''' Calculate no fiber collision weights fraction 

    '''
    # catalog and correction dictionary 
    catalog = cat_corr['catalog']
    correction = cat_corr['correction'] 
    
    if correction['name'] != 'upweight': 
        raise NameError('asdlkfjadf') 

    # import data 
    data = fc_data.galaxy_data('data', clobber=False, **cat_corr) 
    print data.file_name
    
    if catalog['name'].lower() == 'bigmd':  # Big MD 
        wfc = data.wfc 
        comp = np.array([1.0 for i in range(len(wfc))]) 

    elif catalog['name'].lower() in ('nseries', 'qpm'):  # Nseries 
        wfc = data.wfc
        comp = data.comp 

    elif catalog['name'].lower() in ('lasdamasgeo', 'ldgdownnz'):   # LDG 
        wfc = data.weight
        comp = np.array([1.0 for i in range(len(wfc))]) 

    else: 
        raise NotImplementedError('asdlkfjkajsdf') 
   
    w_tot = np.sum(comp) 
    w_withfiber = np.sum(comp[np.where(wfc > 0)]) 
    print 'total weight', w_tot
    print w_tot - w_withfiber, ' galaxies without fibers'
    print 'w_fc = 0 fraction = ', 1.0 - w_withfiber/w_tot

if __name__=='__main__': 
    '''
    corrections = [{'name': 'bigfc'}]
    for i_mock in np.arange(1, 2): 
        for letter in ['a', 'b', 'c', 'd']: 
            for corr in corrections: 
                cat_corr = {
                        'catalog': {'name': 'lasdamasgeo', 'n_mock': i_mock, 'letter': letter}, 
                        'correction': corr} 
                fc_data.galaxy_data('data', clobber=True, **cat_corr) 
    '''
    build_pk('bigmd1', 1, grid=960, quad=False) 
    build_pk('bigmd2', 1, grid=960, quad=False) 
    #build_pk('ldgdownnz', 10, clobber=True, quad=False) 
    #build_pk('lasdamasgeo', 10, clobber=False, grid=360, quad=False) 
    #build_pk('patchy', 10, clobber=False, grid=960, quad=False) 
    #build_pk('cmass', 1, cosmology='fiducial', quad=False)
    #build_pk('nseries', 84, quad=True)
