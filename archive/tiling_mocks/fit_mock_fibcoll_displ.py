import numpy as np
import pylab as py
from scipy.integrate import simps
from scipy.optimize import curve_fit
from matplotlib import rc
import sys

rc('text', usetex=True)
rc('font', family='serif')

dlos_dir = '/mount/riachuelo1/hahn/data/tiling_mocks/'
dlos_file = 'cmass-boss5003sector-icoll012.disp_los.dat'
dlos_data = np.loadtxt(dlos_dir+dlos_file)

mpc_bin = -1000.0+0.1*np.array(range(20001))

fig2 = py.figure(2)
dump = fig2.add_subplot(111)

fig1 = py.figure(1)
ax1 = fig1.add_subplot(111)

hist_dlos = dump.hist(dlos_data, mpc_bin, label='Line of Sight Displacement Histogram')

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
ax1.set_ylim([0,900])
ax1.set_xlabel('Displacement (Mpc)',fontsize=18)
ax1.set_ylabel('Number of Galaxies',fontsize=18)
ax1.legend(loc='best')

for d in [20, 30, 40]:
    RMSfrac = float(len(dlos_data[(dlos_data<d) & (dlos_data>-d)]))/float(len(dlos_data))*100.0
    caption = r''+str(np.int(RMSfrac))+"$\%$"
    ax1.annotate(caption, (float(d),10), xycoords='data', xytext=(float(d), 300), textcoords='data',
            arrowprops=dict(arrowstyle="fancy", facecolor='black', connectionstyle="angle3,angleA=0,angleB=-90"),
            fontsize=20, horizontalalignment='center', verticalalignment='top')
fig1.text(0.5,0.96, 'CMASS Tiling Mocks 5003 Sector icol012', ha='center', va='center', fontsize=20)
py.show()
