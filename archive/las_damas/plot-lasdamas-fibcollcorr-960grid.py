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
file_suffix = '.grid960.P020000.box3600'

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
ax11.loglog( k, avg_delta-avg_true, 'b', linewidth=3,
        label=r"$\overline{P(k)}_{\Delta}-\overline{P(k)}_{\rm{True}}$")
ax11.set_xlim([10**-3,10**0])
ax11.set_ylim([10**2,10**5.5])
ax11.set_xlabel('k',fontsize=20)
ax11.set_ylabel('P(k)',fontsize=20)
ax11.legend(loc='upper right',prop={'size':14})
ax11.grid(True)

ratio_delta_true    = avg_delta/avg_true

ax12 = fig1.add_subplot(122)
ax12.grid(True)
ax12.set_xscale('log')
ax12.set_xlim([10**-3,10**0])

for i in nrange: 
    for letter in ['a','b','c','d']:
        file_peaknbar  = file_prefix+str(i+100)[1:3]+letter+file_midfix+'peaknbar.sigma6.687.fpeak1.0.fibcoll.dat.peaknbar'+file_suffix
        data_peaknbar   = np.loadtxt(dir+file_peaknbar)
    if i==1:
        sum_peaknbar    = data_peaknbar[:,1]
    else:
        sum_peaknbar    += data_peaknbar[:,1]
avg_peaknbar    = sum_peaknbar/float(n)

ax11.loglog( k, avg_peaknbar, color='r', linewidth=2)

ratio_peaknbar_true = avg_peaknbar/avg_true

ax12.plot( k, ratio_peaknbar_true,'r', 
    label=r"$\overline{P(k)}_{\rm{Peak}+\bar{n}(z);\sigma = 6.687,f_{\rm{peak}}=1.0}/\overline{P(k)}_{\rm{True}}$")
ax12.plot( k, ratio_delta_true, 'g',linewidth=2,
        label=r"$\overline{P(k)}_{\Delta}/\overline{P(k)}_{\rm{True}}$")
#ax12.set_ylim([0.85,1.1])
ax12.legend(loc='upper left',prop={'size':14})
fig1.savefig('/home/users/hahn/figures/boss/fiber_collision/lasdamas-fibcoll-peaknbar-delta-comparison-960grid.png',bbox_inches=0)
py.show()

