import numpy as np
import pylab as py
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import interpolate
import sys

rc('text', usetex=True)
rc('font', family='serif')

dir = '/mount/riachuelo1/hahn/power/manera_mock/v11p0/'

file_prefix = 'power_cmass_dr11_north_ir4'
file_midfix = '.v11.0.600randoms.'
file_suffix = '.grid360.P020000.box3600'

n = int(sys.argv[1])
nrange = range(1,n+1)
for i in nrange:
    file_true       = file_prefix+str(i+1000)[1:4]+file_midfix+'wboss'+file_suffix
    file_delta      = file_prefix+str(i+1000)[1:4]+file_midfix+'upweight'+file_suffix
    data_true       = np.loadtxt(dir+file_true)
    data_delta      = np.loadtxt(dir+file_delta)
    if (i==1):
        k               = data_true[:,0]
        sum_true        = data_true[:,1]
        sum_delta       = data_delta[:,1]
    else:
        sum_true        += data_true[:,1]
        sum_delta       += data_delta[:,1]
avg_true        = sum_true/float(n)
avg_delta       = sum_delta/float(n)

fig1 = plt.figure(1, figsize=(7,8))
#ax11 = fig1.add_subplot(121)
#ax11.text(2.0*10**-3.0,10**5.2,str(n)+' PTHalo v11.0 Mocks',fontsize=15)
#ax11.scatter( k, avg_true, color='k',
#        label=r"$\overline{P(k)}_{w_{\rm{BOSS}}}$")
#ax11.loglog( k, avg_delta, 'b', linewidth=3,
#        label=r"$\overline{P(k)}_{\Delta}$")
#ax11.set_xlim([10**-3,10**0])
#ax11.set_ylim([10**3,10**5.5])
#ax11.set_xlabel('k',fontsize=20)
#ax11.set_ylabel('P(k)',fontsize=20)
#ax11.legend(loc='lower left',prop={'size':14})
#ax11.grid(True)

ratio_delta_true    = avg_delta/avg_true

outputname='pthalo-v11p0-upweight-over-wbossonly-ratio.dat'
outputfile=open(outputname,'w')
for j in range(len(k)): 
   outputfile.write(str(k[j])+'\t'+str(ratio_delta_true[j])+'\n') 

ax12 = fig1.add_subplot(111)
ax12.grid(True)
ax12.set_xscale('log')
ax12.set_xlim([10**-3,10**0])

peaknbardir = '/mount/riachuelo1/hahn/power/manera_mock/v11p0/fibcoll/'
for i in nrange: 
    file_peaknbar  = file_prefix+str(i+1000)[1:4]+file_midfix+'peaknbar.sigma5.44fpeak1.0'+file_suffix
    data_peaknbar   = np.loadtxt(peaknbardir+file_peaknbar)
    if (i==1):
        sum_peaknbar    = data_peaknbar[:,1]
    else:
        sum_peaknbar    += data_peaknbar[:,1]
avg_peaknbar    = sum_peaknbar/float(n)

#ax11.loglog( k, avg_peaknbar, color='r', linewidth=2)

ratio_peaknbar_true = avg_peaknbar/avg_true

ax12.text(2.0*10**-3.0,1.2,str(n)+' PTHalo v11.0 Mocks',fontsize=15)
ax12.plot( k, ratio_delta_true, 'g',linewidth=2,
        label=r"$\overline{P(k)}_{\Delta}/\overline{P(k)}_{w_{\rm{BOSS}}}$")
#ax12.plot( k, ratio_peaknbar_true,'r--', linewidth=3, 
#    label=r"$\overline{P(k)}_{\rm{Peak}+\bar{n}(z);\sigma = 5.44,f_{\rm{peak}}=1.0}/\overline{P(k)}_{w_{\rm{BOSS}}}$")
#ax12.set_ylim([0.85,1.1])
ax12.legend(loc='upper left',prop={'size':14})
fig1.savefig('/home/users/hahn/figures/boss/fiber_collision/pthalo-v11p0-fibcoll-delta-wbossonly-comparison-zlim-360grid.png',bbox_inches=0)
py.show()

