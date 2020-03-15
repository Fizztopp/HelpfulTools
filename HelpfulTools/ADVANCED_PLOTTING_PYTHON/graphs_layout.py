# https://matplotlib.org/index.html here you find everything you need

import matplotlib.pyplot as plt  
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 12
mpl.rcParams['font.size'] = 16  # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 16
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.figsize'] = [10.,7.5]

#define some nice color (though most colors are defined in teh standard library)
RED = '#e41a1c'
BLUE = '#377eb8'
GREEN = '#4daf4a'
BROWN = '#fdae61'
VIOLETT = '#6a3d9a' 

# load x data
file_data_x = open('x_data.txt','r')
X = np.loadtxt(file_data_x)
file_data_x.close()

# load y data
file_data_y = open('y_data.txt','r')
Y = np.loadtxt(file_data_y)
file_data_y.close()

# define a useful function 
def gauss_dist(sigma, shift, x):
    return np.exp(-0.5*((x-shift)/sigma)**2)

x_help = np.linspace(1.,10.,100)

# create figure environment
fig = plt.figure(1)
gs1 = gridspec.GridSpec(2, 1) # <- defines grid for subplots


# 1st suplot
ax1 = fig.add_subplot(gs1[0,0])
ax1.set_title(r'title 1',y=1.05)
ax1.set_xlabel(r'$\mathrm{xlabel1}$ $\mathrm{(cm)}$')
ax1.set_ylabel('$\mathrm{ylabel1}$')
ax1.plot(X, Y , 'x', color=BLUE, mew=2, label=r'plotlabel1')
ax1.plot(X, Y , color=RED, label=r'plotlabel2')
ax1.legend(bbox_to_anchor=(0.99, 0.05), loc=4, borderaxespad=0.) 
plt.grid(True)

# 2nd (more advanced) suplot
ax2 = fig.add_subplot(gs1[1,0])
ax2.set_title(r'title 2',y=1.05)
ax2.set_xticks([0, 2, 4, 6, 8, 10])
ax2.set_xticklabels(['$\mathrm{\Phi}$','$\mathrm{\Delta}$','$\mathrm{\Gamma}$','$\mathrm{\Delta}$','$\mathrm{\Xi}$','$\mathrm{\Phi}$'])
ax2.set_xlabel('$\mathrm{xlabel2}$  $\mathrm{(\AA)}$')
ax2.set_ylabel('$\mathrm{ylabel2}$')
ax2.plot(X, Y**2, 'ko', label=r'plotlabel1')
ax2.plot(X, Y**2 , color=VIOLETT, linestyle="--", label=r'plotlabel2')
ax2.fill(x_help,np.amax(Y**2)*gauss_dist(1, 5, x_help), color="gainsboro", facecolor='gainsboro',label=r"$\mathrm{pump}$")
ax2.legend(bbox_to_anchor=(0.95, 0.95), loc=1, borderaxespad=0.) 
ax2.set_xlim(0,10)
plt.grid(True)

#inset for 2nd subplot
ax21 = plt.axes([0,0,1,1])
ip = InsetPosition(ax2, [0.08,0.75,0.2,0.4])

ax21.set_axes_locator(ip)
mark_inset(ax2, ax21, loc1=3, loc2=4, fc="none", ec='0.5')
ax21.plot(X[1:3], Y[1:3]**2,'ko', mew=2, alpha=0.8, label=r'$\mathrm{DFT}$')
ax21.plot(X[1:3], Y[1:3]**2, color=VIOLETT, linestyle='dashed', label=r'$\mathrm{quadratic}$ $\mathrm{fit}$')
ax21.locator_params(nbins=4, axis='x')
ax21.locator_params(nbins=4, axis='y')

plt.tight_layout()
plt.savefig("plot_name" + ".png")
plt.show()


