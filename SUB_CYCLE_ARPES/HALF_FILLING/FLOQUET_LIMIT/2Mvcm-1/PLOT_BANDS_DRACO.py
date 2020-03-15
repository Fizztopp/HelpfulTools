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
mpl.rcParams['figure.figsize'] = [8.,5]

m = 1

RED = '#e41a1c'
BLUE = '#377eb8'
GREEN = '#4daf4a'
BROWN = '#fdae61'
VIOLETT = '#6a3d9a' 


file_BANDS = open('bands0.txt','r')
MAT_BANDS = np.loadtxt(file_BANDS)
file_BANDS.close()

N_BAND = np.size(MAT_BANDS)
print(N_BAND)
k=np.linspace(-np.pi,np.pi,N_BAND)     
       

fig1 = plt.figure(1)
gs1 = gridspec.GridSpec(1, 1)

ax11 = fig1.add_subplot(gs1[0,0])
ax11.set_ylabel(r'$\mathrm{E(k)}$',fontsize=24)
ax11.set_xlabel(r'$\mathrm{ka}$',fontsize=24)
ax11.plot(k,[0]*N_BAND,'k--', linewidth=2.0)
ax11.plot(k,MAT_BANDS, 'k', linewidth=2.0)  
ax11.set_xlim(-np.pi,np.pi)      
plt.tight_layout()
plt.show()
