import numpy as np
import pylab as py
from scipy.integrate import simps
from scipy.optimize import curve_fit
from matplotlib import rc
import sys

rc('text', usetex=True)
rc('font', family='serif')

lasdamas_dir = '/mount/riachuelo1/hahn/data/las_damas/box/'
dlos_file = 'mock_gamma_lrg21p2_zmo_Oriana_1001_z0p342_fof_b0p2.zdist.dlos.dat'
data_file = 'mock_gamma_lrg21p2_zmo_Oriana_1001_z0p342_fof_b0p2.zdist.fibcol.dat' 
dlos = np.loadtxt(lasdamas_dir+dlos_file)
print 'Number of dLOS = ', len(dlos)

mpc_bin = -1000.0+0.1*np.array(range(20001))

fig2 = py.figure(2)
dump = fig2.add_subplot(111)

fig1 = py.figure(1)
ax1 = fig1.add_subplot(111)

hist_dlos = dump.hist(dlos, mpc_bin, label='Line of Sight Displacement Histogram')

def gauss(x,sig):
    return np.max(hist_dlos[0])*np.exp(-0.5*x**2/sig**2)
def expon(x,sig): 
    return np.max(hist_dlos[0])*np.exp(-x/sig)

dlos_x = [ (hist_dlos[1][i] + hist_dlos[1][i+1])/2.0 for i in range(len(hist_dlos[1])-1) ]

popt, pcov = curve_fit(expon, np.array(dlos_x[10000:10500]), hist_dlos[0][10000:10500])
print popt

ax1.plot(dlos_x, hist_dlos[0],linewidth=3, label=r'Histogram of $d_{LOS}$')
ax1.plot(dlos_x[10000:10500], expon(np.array(dlos_x[10000:10500]), popt[0]), 'r', linewidth=3, label=r'Exponential distribution with $\sigma=$'+str(popt))
ax1.set_xlim([-50,50])
ax1.set_xlabel('Displacement (Mpc)')
ax1.set_ylabel('Number of Galaxies')
ax1.legend(loc='best')

for d in [20, 30, 40]:
    RMSfrac = float(len(dlos[(dlos<d) & (dlos>-d)]))/float(len(dlos))*100.0
    caption = r''+str(np.int(RMSfrac))+"$\%$"
    ax1.annotate(caption, (float(d),5), xycoords='data', xytext=(float(d), 200), textcoords='data',
            arrowprops=dict(arrowstyle="fancy", facecolor='black', connectionstyle="angle3,angleA=0,angleB=-90"),
            fontsize=20, horizontalalignment='center', verticalalignment='top')
py.show()
