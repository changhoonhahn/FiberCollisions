'''

Code used to generate results for Sergio's paper 

Author(s) : ChangHoon Hahn 


'''
import numpy as np 
import os.path
import time 
import subprocess
import cosmolopy as cosmos

import fibercollisions as fc
import fibcol_data as fc_data
import plot_fibcol as fc_plot 


def build_pk_comp(option): 
    ''' Calculate P(k) for the following comparison: 

    option == 1:  
        BigMD Nearest Neighbor Upweights with corrected shotnoise
        CMASS Nearest Neighbor Upweights with corrected shotnoise 
    option == 2: 
        BigMD True with corrected shotnoise
        CMASS Fiber Collision corrected

    '''
    if option == 1: 
        kwargs = {'grid': 960, 'quad': False, 'clobber': True}
        fc.build_pk(['bigmd3', 'upweight', 1, kwargs]) 
        fc.build_pk(['cmass', 'upweight', 1, kwargs]) 

def plot_pk_comp(option): 
    if option == 1: 
        catcorr_methods = [
                {'catalog': {'name': 'cmass', 'cosmology': 'fiducial'}, 
                    'correction': {'name': 'upweight'}},
                {'catalog': {'name': 'bigmd3'}, 'correction': {'name': 'upweight'}}
                ]
    elif option == 2: 
        catcorr_methods = [
                {'catalog': {'name': 'cmass', 'cosmology': 'fiducial'}, 
                    'correction': {'name': 'peakshot', 'sigma':6.9, 'fpeak':0.7, 'fit':'gauss'}},
                {'catalog': {'name': 'bigmd3'}, 'correction': {'name': 'true'}}
                ]
    n_mock_list = [1,1]
    fc_plot.plot_pk_fibcol_comp(catcorr_methods, n_mock_list, 
            quad=False, Ngrid=960, type='regular', 
            xrange=[0.001, 1.0], yrange=[10**2, 3*10**5])
    fc_plot.plot_pk_fibcol_comp(catcorr_methods, n_mock_list, 
            quad=False, Ngrid=960, type='ratio', 
            xrange=[0.001, 1.0], yrange=[0.0, 2.0])

if __name__=="__main__": 
    plot_pk_comp(1)
