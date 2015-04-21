import numpy as np
import pylab as py
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import interpolate
import sys

rc('text', usetex=True)
rc('font', family='serif')

dir = '/mount/riachuelo1/hahn/power/las_damas/fiber_collision/'

file_prefix = 'power_sdssmock_gamma_lrgFull_zm_oriana'
file_midfix = '_no.rdcz.'
file_suffix = '.grid360.P020000.box3600'

n = int(sys.argv[1])
nrange = range(1,n+1)
for i in nrange:
    for letter in ['a','b','c','d']:
        file_true       = file_prefix+str(i+100)[1:3]+letter+file_midfix+'dat'+file_suffix
        file_delta      = file_prefix+str(i+100)[1:3]+letter+file_midfix+'fibcoll.dat.delta'+file_suffix
        data_true       = np.loadtxt(dir+file_true)
        data_delta      = np.loadtxt(dir+file_delta)
    if i==1:
        k               = data_true[:,0]
        sum_true        = data_true[:,1]
        sum_delta       = data_delta[:,1]
    else:
        sum_true        += data_true[:,1]
        sum_delta       += data_delta[:,1]
avg_true        = sum_true/float(n)
avg_delta       = sum_delta/float(n)

fig1 = plt.figure(1, figsize=(14,8))
ax11 = fig1.add_subplot(121)
ax11.text(2.0*10**-3.0,10**5.2,str(n)+' Las Damas Mocks',fontsize=15)
ax11.scatter( k, avg_true, color='k',
        label=r"$\overline{P(k)}_{\rm{True}}$")
ax11.loglog( k, avg_delta, 'b', linewidth=3,
        label=r"$\overline{P(k)}_{\Delta}$")
ax11.set_xlim([10**-3,10**0])
ax11.set_ylim([10**3,10**5.5])
ax11.set_xlabel('k',fontsize=20)
ax11.set_ylabel('P(k)',fontsize=20)
ax11.legend(loc='lower left',prop={'size':14})
ax11.grid(True)

ratio_delta_true    = avg_delta/avg_true

ax12 = fig1.add_subplot(122)
ax12.grid(True)
ax12.set_xscale('log')
ax12.set_xlim([10**-3,10**0])

sigma_color = ['#300000','#400000','#500000','#600000','#700000','#800000','#FF0000','#A00000']
for sig in ['0','1','2','3','4','5','6','7']:
    for i in nrange: 
        for letter in ['a','b','c','d']:
            file_peaknbar  = file_prefix+str(i+100)[1:3]+letter+file_midfix+'peaknbar.sigma'+sig+'.687.fpeak1.0.fibcoll.dat.peaknbar'+file_suffix
            data_peaknbar   = np.loadtxt(dir+file_peaknbar)
        if i==1:
            sum_peaknbar    = data_peaknbar[:,1]
        else:
            sum_peaknbar    += data_peaknbar[:,1]
    avg_peaknbar    = sum_peaknbar/float(n)

    ax11.loglog( k, avg_peaknbar, color=sigma_color[int(sig)], linewidth=2)

    ratio_peaknbar_true = avg_peaknbar/avg_true
    if sig=='6': 
        ax12.plot( k, ratio_peaknbar_true, sigma_color[int(sig)], 
            label=r"$\overline{P(k)}_{\rm{Peak}+\bar{n}(z);\sigma = 6.687,f_{\rm{peak}}=1.0}/\overline{P(k)}_{\rm{True}}$")
    else: 
        ax12.scatter( k, ratio_peaknbar_true, color=sigma_color[int(sig)])
    ax12.text(3.2*10**-1.0,ratio_peaknbar_true[-1],r"$f_{\rm{peak}}=1.0, \sigma=$"+sig+"$.687$",fontsize=12)

fpeak_color = ['#3000FF','#4000FF','#5000FF','#6000FF','#7000FF','#8000FF','#9000FF','#A000FF','#B000FF','#C000FF']
for sig in ['0','1','2','3','4','5','6','7','8','9']:
    for i in nrange: 
        for letter in ['a','b','c','d']:
            file_peaknbar  = file_prefix+str(i+100)[1:3]+letter+file_midfix+'peaknbar.sigma6.687.fpeak0.'+sig+'.fibcoll.dat.peaknbar'+file_suffix
            data_peaknbar   = np.loadtxt(dir+file_peaknbar)
        if i==1:
            sum_peaknbar    = data_peaknbar[:,1]
        else:
            sum_peaknbar    += data_peaknbar[:,1]
    avg_peaknbar    = sum_peaknbar/float(n)

    ratio_peaknbar_true = avg_peaknbar/avg_true
    ax12.plot( k, ratio_peaknbar_true, fpeak_color[int(sig)], linewidth=2)
    ax12.text(3.2*10**-1.0,ratio_peaknbar_true[-1],r"$f_{\rm{peak}}=0.$"+sig+'$, \sigma=6.687$',fontsize=12)

ax12.plot( k, ratio_delta_true, 'g',linewidth=2,
        label=r"$\overline{P(k)}_{\Delta}/\overline{P(k)}_{\rm{True}}$")
ax12.set_ylim([0.85,1.1])
ax12.legend(loc='upper left',prop={'size':14})
#fig1.savefig('lasdamas-fibcoll-fpeak_sigma_pk_comparison.png',bbox_inches=0)
py.show()
