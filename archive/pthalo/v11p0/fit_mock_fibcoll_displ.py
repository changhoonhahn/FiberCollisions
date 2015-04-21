import numpy as np
import pylab as py
from scipy.integrate import simps
from scipy.optimize import curve_fit
from matplotlib import rc
import sys

rc('text', usetex=True)
rc('font', family='serif')

n = int(sys.argv[1])
lasdamas_dir = '/mount/riachuelo1/hahn/data/manera_mock/v11p0/fibcoll/'
disp_los = np.loadtxt(lasdamas_dir+'cmass_dr11_north_ir4'+str(n+1000)[1:4]+'.v11.0.wghtv.disp_los.dat')

mpc_bin = -1000.0+0.1*np.array(range(20001))

fig2 = py.figure(2)
dump = fig2.add_subplot(111)

fig1 = py.figure(1)
ax1 = fig1.add_subplot(111)

hist_disp_los = dump.hist(disp_los, mpc_bin, label='Line of Sight Displacement Histogram')

def gauss(x,sig):
    return np.max(hist_disp_los[0])*np.exp(-0.5*x**2/sig**2)
def expon(x,sig): 
    return np.max(hist_disp_los[0])*np.exp(-x/sig)

disp_los_x = [ (hist_disp_los[1][i] + hist_disp_los[1][i+1])/2.0 for i in range(len(hist_disp_los[1])-1) ]

popt, pcov = curve_fit(expon, np.array(disp_los_x[10000:10500]), hist_disp_los[0][10000:10500])
print popt

ax1.plot(disp_los_x, hist_disp_los[0],linewidth=3, label=r'Histogram of $d_{LOS}$')
ax1.plot(disp_los_x[10000:10500], expon(np.array(disp_los_x[10000:10500]), popt[0]), 'r', linewidth=3, label=r'Exponential distribution with $\sigma=$'+str(popt))
ax1.set_xlim([-50,50])
ax1.set_ylim([0,200])
ax1.set_xlabel('Displacement (Mpc)',fontsize=18)
ax1.set_ylabel('Number of Galaxies',fontsize=18)
ax1.legend(loc='best')

for d in [20, 30, 40]:
    RMSfrac = float(len(disp_los[(disp_los<d) & (disp_los>-d)]))/float(len(disp_los))*100.0
    caption = r''+str(np.int(RMSfrac))+"$\%$"
    ax1.annotate(caption, (float(d),5), xycoords='data', xytext=(float(d), 25), textcoords='data',
            arrowprops=dict(arrowstyle="fancy", facecolor='black', connectionstyle="angle3,angleA=0,angleB=-90"),
            fontsize=20, horizontalalignment='center', verticalalignment='top')
fig1.text(0.5,0.96, 'PTHalo v11.0 ir4'+str(n+1000)[1:4]+' North', ha='center', va='center', fontsize=20)
py.show()
