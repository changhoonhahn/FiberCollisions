'''

Convolving the window function to the Nseries Box power spectra

'''

def W_k_kp_2(k, kp, l, L): 
    '''

    |W(k, k')|^2_{l, L}

    '''
    alpha = 2. * (1j)**l * (-1j)**L * (2.*l+1)
    alpha_real = alpha.real

    # read in random 

    #

