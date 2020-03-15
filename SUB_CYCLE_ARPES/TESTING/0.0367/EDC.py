import matplotlib.pyplot as plt  
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 15
mpl.rcParams['font.size'] = 20  # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 24
mpl.rcParams['axes.titlesize'] = 24
mpl.rcParams['axes.labelsize'] = 24
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.titlesize'] = 28
mpl.rcParams['figure.figsize'] = [10,10]
mpl.rcParams['image.cmap'] = 'seismic'

tc = 0.658 # 1/eV = 0.658 fs

RED = '#e41a1c'
BLUE = '#377eb8'
GREEN = '#4daf4a'
BROWN = '#fdae61'
VIOLETT = '#6a3d9a' 

def gauss(sigma, shift, x):
    return np.exp(-0.5*((x-shift)/sigma)**2)

fac = 4/3.84*1/10 # 10 meV / 3.4 AA = 1MV/cm  

AA = 10.0*fac
omega_min = -0.9
omega_max = +0.1
N_OMEGA = 128

OMEGA = np.linspace(omega_min, omega_max, N_OMEGA)

file_EDC = open('DATA/IPHOTO10.txt','r')
EDC_1e1 = np.loadtxt(file_EDC)
file_EDC.close()
file_EDC = open('DATA/IPHOTO100.txt','r')
EDC_1e2 = np.loadtxt(file_EDC)
file_EDC.close()
file_EDC = open('DATA/IPHOTO1000.txt','r')
EDC_1e3 = np.loadtxt(file_EDC)
file_EDC.close()
file_EDC = open('DATA/IPHOTO10000.txt','r')
EDC_1e4 = np.loadtxt(file_EDC)
file_EDC.close()

file_EDC = open('DATA/IPHOTO6.453320.txt','r')
EDC_10 = np.loadtxt(file_EDC)
file_EDC.close()
file_EDC = open('DATA/IPHOTO12.906600.txt','r')
EDC_20 = np.loadtxt(file_EDC)
file_EDC.close()
file_EDC = open('DATA/IPHOTO25.813300.txt','r')
EDC_40 = np.loadtxt(file_EDC)
file_EDC.close()
file_EDC = open('DATA/IPHOTO38.719900.txt','r')
EDC_60 = np.loadtxt(file_EDC)
file_EDC.close()
file_EDC = open('DATA/IPHOTO51.626600.txt','r')
EDC_80 = np.loadtxt(file_EDC)
file_EDC.close()
file_EDC = open('DATA/IPHOTO64.533200.txt','r')
EDC_100 = np.loadtxt(file_EDC)
file_EDC.close()


fig1 = plt.figure(1)
gs1 = gridspec.GridSpec(1, 1)

ax00 = fig1.add_subplot(gs1[0,0])
ax00.set_title("EDC @$k=0.0$, $\mathrm{\Delta t}=0.5$ $\mathrm{fs}$")
ax00.set_xlim(omega_min, omega_max)
ax00.set_xticks([-0.8, -0.6, -0.4, -0.2, -0.0])
#ax01.set_yticklabels(["-1","0","+1"])

#ax01.set_ylim(-1., 1.)
ax00.set_ylabel('$\mathrm{I_{photo}}$')
ax00.set_xlabel('$\mathrm{Energy}$ $(\mathrm{eV})$')
ax00.plot(OMEGA, EDC_10, color='k', label=r"$\mathrm{FWHM}=10$ $\mathrm{fs}$")
ax00.plot(OMEGA, EDC_20, color=BROWN, label=r"$\mathrm{FWHM}=20$ $\mathrm{fs}$")
ax00.plot(OMEGA, EDC_40, color=RED, label=r"$\mathrm{FWHM}=40$ $\mathrm{fs}$")
ax00.plot(OMEGA, EDC_60, color=BLUE, label=r"$\mathrm{FWHM}=60$ $\mathrm{fs}$")
ax00.plot(OMEGA, EDC_80, color=VIOLETT, label=r"$\mathrm{FWHM}=80$ $\mathrm{fs}$")
ax00.plot(OMEGA, EDC_100, color=GREEN, label=r"$\mathrm{FWHM}=100$ $\mathrm{fs}$")
plt.legend()


fig1 = plt.figure(2)
gs1 = gridspec.GridSpec(1, 1)

ax01 = fig1.add_subplot(gs1[0,0])
ax01.set_title("EDC @$k=0.0$, $\mathrm{FWHM}=100$ $\mathrm{fs}$")
ax01.set_xlim(omega_min, omega_max)
ax01.set_xticks([-0.8, -0.6, -0.4, -0.2, -0.0])
#ax01.set_yticklabels(["-1","0","+1"])

#ax01.set_ylim(-1., 1.)
ax01.set_ylabel('$\mathrm{I_{photo}}$')
ax01.set_xlabel('$\mathrm{Energy}$ $(\mathrm{eV})$')
ax01.plot(OMEGA, EDC_1e1, 'k--', linewidth=0.5, label=r"$\mathrm{\Delta t}=50$ $\mathrm{fs}$")
ax01.plot(OMEGA, EDC_1e2, color=RED, label=r"$\mathrm{\Delta t}=5$ $\mathrm{fs}$")
ax01.plot(OMEGA, EDC_1e3, color=BLUE, label=r"$\mathrm{\Delta t}=0.5$ $\mathrm{fs}$")
ax01.plot(OMEGA, EDC_1e4, color=VIOLETT, label=r"$\mathrm{\Delta t}=0.05$ $\mathrm{fs}$")
plt.legend()

plt.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.5)
plt.show()  
