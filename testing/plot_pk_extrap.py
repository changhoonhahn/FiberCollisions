'''
plotting for pk_extrap

'''
import pickle
import numpy as np

import pk_extrap 
from defutility.plotting import prettyplot 
from defutility.plotting import prettycolors 

def plot_avg_Plk(l, n_mocks, Ngrid=360, yscale='linear', **kwargs): 
    '''
    Plot average P_l(k)
    '''
    if not isinstance(l, list): 
        ls = [l]
    else: 
        ls = l 

    prettyplot()
    pretty_colors = prettycolors()

    # plot P(k) data
    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    for l_i in ls: 
        k_i, plk_i = pk_extrap.average_Pk(l_i, n_mocks, Ngrid=Ngrid)
        
        sub.plot(k_i, plk_i, label=r"$l = "+str(l_i)+"$", c=pretty_colors[l_i+1], lw=4)
    
    if 'xrange' in kwargs.keys(): 
        sub.set_xlim(kwargs['xrange'])
    else:
        sub.set_xlim([1e-3, 1.0])

    if 'yrange' in kwargs.keys(): 
        sub.set_ylim(kwargs['yrange'])
    else: 
        pass

    sub.set_ylabel(r"$\mathtt{P_l(k)}$", fontsize=40)
    sub.set_xlabel(r"$\mathtt{k}$", fontsize=40)

    sub.set_xscale('log')
    if yscale == 'log': 
        sub.set_yscale('log')
    else: 
        pass 

    sub.legend(loc='upper right')

    fig_file = ''.join([
        'qaplot_avgP', ''.join([str(l_i) for l_i in l]), 'k_',
        'Ngrid', str(Ngrid), '.png'
        ])
    fig.savefig(fig_file, bbox_inches='tight')
    #plt.show()

def plot_avg_P4k_comp(n_mocks, Ngrid=360, yscale='log', **kwargs): 
    '''
    Plot average P_l(k)
    '''

    prettyplot()
    pretty_colors = prettycolors()

    # plot P(k) data
    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    for i_l, l_i in enumerate([4, '4fast']):
        k_i, plk_i = pk_extrap.average_Pk(l_i, n_mocks, Ngrid=Ngrid)

        l_i_str =  str(l_i)
        try:
            l_i_str =  str(l_i)
        except TypeError: 
            l_i_str = l_i

        sub.plot(
                k_i, 
                np.abs(plk_i), 
                label=r"$l = "+l_i_str+"$", 
                c=pretty_colors[i_l*2+1], 
                lw=4
                )
    
    if 'xrange' in kwargs.keys(): 
        sub.set_xlim(kwargs['xrange'])
    else:
        sub.set_xlim([1e-3, 1.0])

    if 'yrange' in kwargs.keys(): 
        sub.set_ylim(kwargs['yrange'])
    else: 
        pass

    sub.set_ylabel(r"$\mathtt{|P_4(k)|}$", fontsize=40)
    sub.set_xlabel(r"$\mathtt{k}$", fontsize=40)

    sub.set_xscale('log')
    if yscale == 'log': 
        sub.set_yscale('log')
    else: 
        pass 

    sub.legend(loc='upper right')

    fig_file = ''.join([
        'qaplot_avgP4_4fast_k_comparison_',
        'Ngrid', str(Ngrid), '.png'
        ])
    fig.savefig(fig_file, bbox_inches='tight')

def plot_avg_Plk_negative(l, n_mocks, Ngrid=360): 
    '''
    Plot the negative portion of the average P_l(k). The y-axis is linear not logarithmic 
    '''

    if not isinstance(l, list): 
        ls = [l]
    else: 
        ls = l 

    prettyplot()
    pretty_colors = prettycolors()

    # plot P(k) data
    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    for l_i in ls: 
        k_i, plk_i = pk_extrap.average_Pk(l_i, n_mocks, Ngrid=Ngrid)

        negative = np.where(plk_i < 0.)
        
        if l_i < 4: 
            sub.plot(
                    k_i[negative], 
                    plk_i[negative], 
                    label=r"$P_"+str(l_i)+"(k)$", c=pretty_colors[l_i+1], lw=4)
        else:
            sub.scatter(
                    k_i[negative], 
                    plk_i[negative], 
                    label=r"$P_"+str(l_i)+"(k)$", c=pretty_colors[l_i+1])

    sub.set_xlim([1e-3, 1.0])
    sub.set_ylabel(r"$\mathtt{P_l(k)}$", fontsize=40)
    sub.set_xlabel(r"$\mathtt{k}$", fontsize=40)
    sub.set_xscale('log')
    sub.legend(loc='upper right')
    plt.show()

def plot_avg_Plk_extrap(l, n_mocks, Ngrid=360, k_fixed=0.6, k_max=0.3, **kwargs): 
    '''
    Plot average P_l(k) along with extrapolated P_l(k) for k > k_max
    '''
    if not isinstance(l, list): 
        ls = [l]
    else: 
        ls = l 

    prettyplot()
    pretty_colors = prettycolors()

    # plot P(k) data
    fig = plt.figure(figsize=(10,10))
    sub = fig.add_subplot(111)

    for l_i in ls: 

        # average P_l(k)
        k_i, plk_i = pk_extrap.average_Pk(l_i, n_mocks, Ngrid=Ngrid)
        sub.plot(k_i, np.abs(plk_i), label=r"$l = "+str(l_i)+"$", c=pretty_colors[l_i+1], lw=4)

        # extrapolated P_l(k)
        if isinstance(k_max, list) or isinstance(k_max, np.ndarray): 
            for k_max_i in k_max: 
                bestfit_param = pk_extrap.pk_bestfit(k_i, plk_i, k_max=k_max_i, k_fixed=k_fixed)
                if k_max_i == k_max[0]:
                    plot_label = r"$l = "+str(l_i)+"$ extrap."
                else: 
                    plot_label = None

                sub.plot(
                        np.arange(k_max_i, 1.0, 0.05), 
                        np.abs(
                            pk_extrap.pk_powerlaw(
                                np.arange(k_max_i, 1.0, 0.05), bestfit_param, k_fixed=k_fixed
                                )
                            ),
                        label=plot_label, 
                        c='k', #c=pretty_colors[l_i+1], 
                        ls='--', 
                        lw=2
                        )
        else:
            bestfit_param = pk_extrap.pk_bestfit(k_i, plk_i, k_max=k_max, k_fixed=k_fixed)
            print bestfit_param
            sub.plot(
                    np.arange(k_max, 1.0, 0.05), 
                    pk_extrap.pk_powerlaw(
                        np.arange(k_max, 1.0, 0.05), bestfit_param, k_fixed=k_fixed
                        ),
                    label=r"$l = "+str(l_i)+"$ extrap.", 
                    c=pretty_colors[l_i+1], 
                    ls='--', 
                    lw=2
                    )
    
    if 'xrange' in kwargs.keys(): 
        sub.set_xlim(kwargs['xrange'])
    else:
        sub.set_xlim([1e-3, 1.0])

    if 'yrange' in kwargs.keys(): 
        sub.set_ylim(kwargs['yrange'])
    else: 
        pass

    sub.set_ylabel(r"$\mathtt{log\;|P_l(k)|}$", fontsize=40)
    sub.set_xlabel(r"$\mathtt{k}$", fontsize=40)

    sub.set_xscale('log')
    sub.set_yscale('log')

    sub.legend(loc='upper right')

    fig_file = ''.join([
        'qaplot_avgP', ''.join([str(l_i) for l_i in l]), 'k_',
        'Ngrid', str(Ngrid), '_extrapolation.png'
        ])
    fig.savefig(fig_file, bbox_inches='tight')
    #plt.show()

def plot_avg_Pkfast_test(**kwargs): 
    '''
    '''

    prettyplot()
    pretty_colors = prettycolors()


    for i_l, l_i in enumerate([0, 2, 4]):
        # P(k) fast
        k_fast_i, plk_fast_i = pk_extrap.average_Pk(l_i, 1, Ngrid=960)
        k_i, plk_i = pk_extrap.average_Pk(l_i, 1, Ngrid=720)

        l_i_str =  str(l_i)
        try:
            l_i_str =  str(l_i)
        except TypeError: 
            l_i_str = l_i
    
        # plot P(k) data
        fig = plt.figure(figsize=(10,10))
        sub = fig.add_subplot(111)

        sub.scatter(
                k_fast_i, 
                np.abs(plk_fast_i), 
                c=pretty_colors[i_l*2+1]
                )

        sub.plot(
                k_i, 
                np.abs(plk_i), 
                label=r"$l = "+l_i_str+"$", 
                c='k',
                lw=2, 
                ls='--'
                )
    
        if 'xrange' in kwargs.keys(): 
            sub.set_xlim(kwargs['xrange'])
        else:
            sub.set_xlim([1e-3, 1.0])

        sub.set_ylabel(r"$\mathtt{|P_l(k)|}$", fontsize=40)
        sub.set_xlabel(r"$\mathtt{k}$", fontsize=40)

        sub.set_xscale('log')
        sub.set_yscale('log')

        sub.legend(loc='upper right')

        fig_file = ''.join(['qaplot_avgP'+str(l_i)+'k_fast_test.png'])
        fig.savefig(fig_file, bbox_inches='tight')
        plt.close()

if __name__=="__main__":
    plot_avg_Pkfast_test(xrange=[0.1, 1.0])
    #plot_avg_P4k_comp(10, Ngrid=720)
    
    #for l_i in [0, 2, 4]:
    #    plot_avg_Plk_extrap(
    #            [l_i], 
    #            10, 
    #            Ngrid=720, 
    #            k_fixed=0.6, 
    #            k_max=np.arange(0.4, 0.65, 0.05),
    #            xrange=[0.1, 1.0], yrange=[10.**0., 10.**5.]
    #            )
