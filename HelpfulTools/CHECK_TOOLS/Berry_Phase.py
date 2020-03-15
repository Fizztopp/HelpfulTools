# Berry phase for effective Monolayer Graphene (for details see: A short Course on Topologiacal Insulators, Sec. 2.5.3)
"""
Created on Mon Mar  4 09:59:22 2019

@author: toppgabr
"""

import numpy as np 
from scipy.linalg import norm         
import matplotlib as mpl
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 15
mpl.rcParams['font.size'] = 20  # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] =15
mpl.rcParams['axes.titlesize'] = 24
mpl.rcParams['axes.labelsize'] = 22
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.titlesize'] = 24
mpl.rcParams['figure.figsize'] = [6.,6]

shift = 0.0
mass = 0.05

s0 = np.eye(2)
s1 = np.array([[0., 1.],[1., 0.]])
s2 = np.array([[0.,-1j],[1j, 0.]])
s3 = np.array([[1., 0.],[0.,-1.]])

def Hk(kx, ky):
   #return shift*s0 + np.cos(kx)*s1 - np.sin(ky)*s2 + mass*s3
    return shift*s0 + kx*s1 + ky*s2 + mass*s3

def Hkdkx(kx, ky):
    #return -np.sin(kx)*s1
    return s1
    
def Hkdky(kx, ky):
    #return -np.cos(kx)*s2
    return s2

def diag(kx, ky):
    eigenvals, eigenvecs = np.linalg.eigh(Hk(kx, ky))
    return eigenvals, eigenvecs ## columns are eigenvectors

def Berry_Curv(kx, ky, n):
    BCz = 0.0
    evals, S = diag(kx, ky)    
    Hdx = Hkdkx(kx, ky)
    Hdy = Hkdky(kx, ky)
    for m in range(2):
        if m==n:
            continue      
        BCz = BCz + (np.dot(np.conj(S[:,n]),np.dot(Hdx,S[:,m]))*np.dot(np.conj(S[:,m]),np.dot(Hdy,S[:,n])) - np.dot(np.conj(S[:,n]),np.dot(Hdy,S[:,m]))*np.dot(np.conj(S[:,m]),np.dot(Hdx,S[:,n])))/(evals[n]-evals[m])**2 
    return -np.imag(BCz)    
#BCz = BCz + (np.dot(S[n,:],np.dot(Hdx,SCT[:,m]))*np.dot(S[m,:],np.dot(Hdy,SCT[:,n])) - np.dot(S[n,:],np.dot(Hdy,SCT[:,m]))*np.dot(S[m,:],np.dot(Hdx,SCT[:,n])))/(evals[n]-evals[m])**2 

## Use Berry's formula
print(Berry_Curv(0.0, 0.0, 0))
print(Berry_Curv(0.0, 0.0, 1))

FAC = 1.0
NN = 100
kx = np.linspace(-FAC*np.pi, +FAC*np.pi, NN) 
ky = np.linspace(-FAC*np.pi, +FAC*np.pi, NN) 

evals0 = np.zeros(NN)
evals1 = np.zeros(NN)

for i in range(NN):
    evals0[i] = diag(kx[i], 0.0)[0][0]
    evals1[i] = diag(kx[i], 0.0)[0][1]

#==============================================================================
# fig0 = plt.figure(3)
# gs0 = gridspec.GridSpec(1, 1)
# ax00 = fig0.add_subplot(gs0[0,0])
# ax00.set_ylabel(r'kx')
# ax00.set_ylabel(r'energy')
# ax00.plot(kx, evals0, 'r-')
# ax00.plot(kx, evals1, 'b-')
# #ax00.set_xlim(-np.pi, np.pi)   
# plt.tight_layout()
# plt.show()
#==============================================================================
AREA = 0
karray = np.zeros((NN*NN,2)) 
X, Y = np.meshgrid(kx, ky)
Z = np.zeros((NN,NN))
temp=0
for i in range(NN):
    for j in range(NN):
        Z[i,j] = Berry_Curv(kx[i], ky[j], 0)
        AREA += (2*FAC*np.pi/(NN))**2
        temp += Z[i,j]#*(2*FAC*np.pi/(NN))**2
        karray[i+NN*j,:] = np.array([kx[i],ky[j]])  
print('Chern number (dH/dk): '+str(-1/(2*np.pi)*temp*(2*FAC*np.pi/(NN))**2))       
print(AREA)
print((2*FAC*np.pi)**2)
#plt.scatter(karray[:,0], karray[:,1], c="r", marker="o", label="basis 1")
#plt.show()

Chern = -1/(2*np.pi)*np.trapz(np.trapz(Z,dx=2*FAC*np.pi/(NN)),dx=2*FAC*np.pi/(NN))
print('Chern number (check): '+str(Chern))


## Calculate Berry Curvature by Sum of local fluxes
#==============================================================================
# fig1 = plt.figure(4)
# gs1 = gridspec.GridSpec(1, 1)
# ax11 = fig1.add_subplot(gs1[0,0])
# ax11.set_xlabel(r'kx')
# ax11.set_ylabel(r'ky')
# ax11.contour(X, Y, Z, 20, cmap='RdYlBu')
# ax11.set_xlim(-np.pi, np.pi)
# ax11.set_ylim(-np.pi, np.pi)     
# plt.tight_layout()
# plt.show()
#==============================================================================

n = 0
kmin = -1e-6*np.pi 
kmax = +1e-6*np.pi
Nk = 2                                                                      # number of plaquettes not important without Gauge field!
kvec = np.array([0.0,0.0]) 
k0 = np.array([0.,0.])
dk = (kmax-kmin)/(Nk-1);
karray = np.zeros((Nk*Nk,2)) 
S_ARRAY = np.zeros((Nk,Nk,2,2),dtype=complex) 
BC_ARRAY = np.zeros((Nk,Nk)) 

def vec_prod(v1, v2):
    size = np.size(v1)
    result = 0.0
    for i in range(size):
        result += np.conj(v1[i])*v2[i]
    return result   
	
k0[0] = kvec[0]-0.5*(kmax-kmin)
k0[1] = kvec[1]-0.5*(kmax-kmin)

for i in range(Nk):
    kvec[0] = k0[0]+i*dk
    #print(kvec[0])
    for j in range(Nk):
        kvec[1] = k0[1]+j*dk
        evals, S_ARRAY[i,j,:,:] = diag(kvec[0], kvec[1])
        karray[i+Nk*j,:] = np.array([kvec[0],kvec[1]]) 

#plt.scatter(karray[:,0], karray[:,1], c="r", marker="o", label="basis 1")
#plt.show()

Bphase = 0.0
CHECK = 0.0
Bflux = 0.0
for i in range(Nk-1):
    for j in range(Nk-1):
        ## sum_nm F_nm                
        Bphase += np.imag(np.log(vec_prod(S_ARRAY[i,j,:,n], S_ARRAY[i+1,j,:,n])*vec_prod(S_ARRAY[i+1,j,:,n], S_ARRAY[i+1,j+1,:,n])*vec_prod(S_ARRAY[i+1,j+1,:,n], S_ARRAY[i,j+1,:,n])*vec_prod(S_ARRAY[i,j+1,:,n], S_ARRAY[i,j,:,n]))) 
        CHECK += np.imag(np.log(vec_prod(S_ARRAY[i,j,:,n], S_ARRAY[i+1,j,:,n])))+np.imag(np.log(vec_prod(S_ARRAY[i+1,j,:,n], S_ARRAY[i+1,j+1,:,n])))+np.imag(np.log(vec_prod(S_ARRAY[i+1,j+1,:,n], S_ARRAY[i,j+1,:,n])))+np.imag(np.log(vec_prod(S_ARRAY[i,j+1,:,n], S_ARRAY[i,j,:,n])))
        print(vec_prod(S_ARRAY[i,j,:,n], S_ARRAY[i+1,j,:,n]*vec_prod(S_ARRAY[i+1,j,:,n], S_ARRAY[i+1,j+1,:,n])*vec_prod(S_ARRAY[i+1,j+1,:,n], S_ARRAY[i,j+1,:,n])*vec_prod(S_ARRAY[i,j+1,:,n], S_ARRAY[i,j,:,n])))
        
print('Berry Curvature (dn/dk): '+str(Bphase/(kmax-kmin)**2)) # kmin,kmax --> 0.0 gives Berry curvature around 0
print('CKECK: '+str(CHECK/(kmax-kmin)**2)) # kmin,kmax --> 0.0 gives Berry curvature around 0
print('Chern number (dn/dk): '+str(Bphase/(2*np.pi)))      

