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
mpl.rcParams['image.cmap'] = 'nipy_spectral'

tc = 0.658 # 1/eV = 0.658 fs

RED = '#e41a1c'
BLUE = '#377eb8'
GREEN = '#4daf4a'
BROWN = '#fdae61'
VIOLETT = '#6a3d9a' 

def gauss(sigma, shift, x):
    return np.exp(-0.5*((x-shift)/sigma)**2)

fac = 207./3.84*1/10 # 10 meV / 3.84 AA = 1MV/cm  

AA = 0.185*fac
omega_min = -1.0
omega_max = +0.5
k_min = 0
k_max = 128
t_min = 0
t_max = 760.0
TIMESTEPS = 1e3
T_PROBE = 380.0                                           		
SIGMA_PROBE  = 65.0   
mu = 0.0

file_ARPES= open('ARPES.txt','r')
ARPES = np.loadtxt(file_ARPES)
file_ARPES.close()

file_Driving = open('DRIVING_t.txt','r')
Driving = np.loadtxt(file_Driving)*fac
file_Driving.close()

file_BANDS = open('bands0.txt','r')
MAT_BANDS = np.loadtxt(file_BANDS)
file_BANDS.close()

N_BAND = np.size(MAT_BANDS)
k=np.linspace(-np.pi,np.pi,N_BAND) 

x_min = -np.pi
x_max = +np.pi
     
t = np.linspace(t_min, t_max, TIMESTEPS)*tc 

fig1 = plt.figure(4)
gs1 = gridspec.GridSpec(4, 4)

ax01 = fig1.add_subplot(gs1[0,0:4])
ax01.plot(t,Driving, "k", linewidth=2.0, label=r"$\mathrm{\Omega = 0.2}$ $\mathrm{eV}$")

ax01.fill(t,-AA+2*AA*gauss(SIGMA_PROBE*tc, T_PROBE*tc, t), color="gainsboro", facecolor='gainsboro',label=r"$\mathrm{probe}$")
#ax01.axvline(x=T_PROBE*tc, ymin=0.0, ymax = 1.0, c="b", linestyle="-")

ax01.set_xlim(t_min,t_max*tc)
ax01.set_ylim(-1.,+1)
ax01.set_ylabel('$\mathrm{A(t)}$ $\mathrm{MV/cm}$')
ax01.set_xlabel('$\mathrm{time}$ $(\mathrm{fs})$')
ax01.xaxis.tick_top()
ax01.xaxis.set_label_position('top')
#plt.legend()


ax00 = fig1.add_subplot(gs1[1:4,0:4])
#ax00.set_title(r'$\mathrm{ARPES}$')
 
#for ii in range(np.size(ARPES[:,0])):
#    for jj in range(np.size(ARPES[0,:])):
#        if(ARPES[ii,jj]>0.2):
#            ARPES[ii,jj]=0.0

MAX = np.amax(ARPES)

HEX=ax00.imshow(ARPES/MAX, aspect='auto', extent=[x_min,x_max ,omega_min, omega_max], vmin=0.0, vmax=1.0)
cbar = fig1.colorbar(HEX, ticks=[0.0, 1.0])
cbar.ax.set_yticklabels(["0","1"])

ax00.set_xticks([-np.pi, 0.0, +np.pi])
ax00.set_xticklabels([r"$-\mathrm{\pi/a}$","$\mathrm{\Gamma}$",r"$+\mathrm{\pi/a}$"], fontsize=24)
ax00.plot(k,[0]*N_BAND,'w--', linewidth=2.0, label=r"$\mathrm{\mu}$")
ax00.plot(k,MAT_BANDS, 'y', linewidth=1.0, label=r"$\mathrm{initial}$ $\mathrm{bands}$")
#plt.legend(loc="lower right")

#ax00.set_xlabel('$\mathrm{momentum}$')
ax00.set_ylabel('$\mathrm{E(k)}$ $(\mathrm{eV})$')
ax00.set_xlim(x_min,x_max)
ax00.set_ylim(omega_min, omega_max)

plt.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.5)
plt.show()  
