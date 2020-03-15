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
mpl.rcParams['image.cmap'] = 'gnuplot'

tc = 0.658 # 1/eV = 0.658 fs

RED = '#e41a1c'
BLUE = '#377eb8'
GREEN = '#4daf4a'
BROWN = '#fdae61'
VIOLETT = '#6a3d9a' 

def gauss(sigma, shift, x):
    return np.exp(-0.5*((x-shift)/sigma)**2)


T_PROBES = np.array([759.87, 949.84, 1139.82, 1329.79, 1519.76, 1709.73, 1899.70, 2089.67, 2279.64, 2469.6, 2659.57, 2849.54, 3039.51, 3229.48, 3419.45, 3609.42, 3799.39])
T_PROBESS = np.array(["759.87", "949.84", "1139.82", "1329.79", "1519.76", "1709.73", "1899.70", "2089.67", "2279.64", "2469.60", "2659.57", "2849.54", "3039.51", "3229.48", "3419.45", "3609.42", "3799.39"])

a = 2
AA = 0.1
omega_min = -0.9
omega_max = +0.1
k_min = 0
k_max = 256
t_min = 0
t_max = 4559.27
TIMESTEPS = 1e3
T_PROBE = 759.87                                       		
SIGMA_PROBE  = 150.0  

for i in range(17):   
    
    file_BANDS = open('mu.txt','r')
    mu  = np.loadtxt(file_BANDS)
    file_BANDS.close()
    
    file_BANDS = open('bands0.txt','r')
    MAT_BANDS = np.loadtxt(file_BANDS)-mu
    file_BANDS.close()
    
    file_ARPES= open('ARPES_'+T_PROBESS[i]+'0000.txt','r')
    ARPES= np.loadtxt(file_ARPES)
    file_ARPES.close()
    
    file_Driving = open('DRIVING_t.txt','r')
    Driving = np.loadtxt(file_Driving)
    file_Driving.close()
    
    
    N_BAND = np.size(MAT_BANDS[:,0])
    k=np.linspace(0,np.pi,N_BAND) 
    
    x_min = 0.0
    x_max = +np.pi
         
    t = np.linspace(t_min, t_max, TIMESTEPS)*tc 
    
    fig1 = plt.figure(2)
    gs1 = gridspec.GridSpec(4, 4)
    
    
    ax01 = fig1.add_subplot(gs1[0,0:4])
    ax01.plot(t,AA*Driving[:,0]/np.amax(Driving[:,0]), "k", linewidth=2.0)
    
    ax01.fill(t,2*AA*gauss(SIGMA_PROBE*tc, T_PROBES[i]*tc, t)-AA, color="gainsboro", facecolor='gainsboro')
    #ax01.axvline(x=T_PROBE*tc, ymin=0.0, ymax = 1.0, c="b", linestyle="-")
    
    ax01.set_xlim(t_min,t_max*tc)
    #ax01.set_ylim(-AA,+AA)
    ax01.set_ylabel('$\mathrm{A_x(t)}$ $(\mathrm{MV cm^{-1}})$')
    ax01.set_xlabel('$\mathrm{time}$ $(\mathrm{fs})$')
    ax01.xaxis.tick_top()
    ax01.xaxis.set_label_position('top')
    plt.legend()
    
    
    ax00 = fig1.add_subplot(gs1[1:4,0:4])
    #ax00.set_title(r'$\mathrm{ARPES}$')
    
    MAX = np.amax(ARPES)
    HEX=ax00.imshow(ARPES/MAX, aspect='auto', extent=[x_min,x_max ,omega_min-mu,omega_max-mu], vmin=0., vmax=1.0)
    plt.colorbar(HEX, ticks=[0.0, 1.0])
    
    ax00.set_xticks([x_min, np.pi/2])
    ax00.set_xticklabels(['$\Gamma$', 'X'], fontsize=24)
    ax00.plot(k,[0]*N_BAND,'w--', linewidth=2.0)
    ax00.plot(k,MAT_BANDS[:,0], 'g', linewidth=2.0, label=r"$\mathrm{initial}$ $\mathrm{bands}$")
    ax00.plot(k,MAT_BANDS[:,1], 'g', linewidth=2.0)
    ax00.plot(k,MAT_BANDS[:,2], 'g', linewidth=2.0)
    ax00.plot(k,MAT_BANDS[:,3], 'g', linewidth=2.0)
    ax00.plot(k,MAT_BANDS[:,4], 'g', linewidth=2.0)
    ax00.plot(k,MAT_BANDS[:,5], 'g', linewidth=2.0)
    ax00.plot(k,MAT_BANDS[:,6], 'g', linewidth=2.0)
    ax00.plot(k,MAT_BANDS[:,7], 'g', linewidth=2.0)
    plt.legend(loc='lower right')
    
    #ax00.set_xlabel('$\mathrm{momentum}$')
    ax00.set_ylabel('$\mathrm{E(k)}$ $(\mathrm{eV})$')
    ax00.set_xlim(x_min,np.pi/2)
    ax00.set_ylim(omega_min-mu,omega_max-mu)

    plt.tight_layout(pad=1.0, w_pad=0.1, h_pad=0.5)
    fig1.savefig('MOVIE/'+str(i)+'.png')
    plt.clf()
plt.close()