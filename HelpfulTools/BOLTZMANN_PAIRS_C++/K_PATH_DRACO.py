# Python programm to calculate the k-points of the full/irreducible Brilluoin zone. Here for the example of an fcc lattice with a four atomic basis and 4 different atom sort
# (c) Gabriel E. Topp, gabriel.topp@mpsd.mpg.de

import numpy as np             
import spglib

mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['lines.linewidth'] = 3
mpl.rcParams['lines.markersize'] = 15
mpl.rcParams['font.size'] = 14  # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 20
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['axes.labelsize'] = 20
mpl.rcParams['xtick.major.size'] = 10
mpl.rcParams['ytick.major.size'] = 10
mpl.rcParams['xtick.major.width'] = 1
mpl.rcParams['ytick.major.width'] = 1
mpl.rcParams['xtick.direction'] = 'inout'
mpl.rcParams['ytick.direction'] = 'inout'
mpl.rcParams['figure.titlesize'] = 28
mpl.rcParams['figure.figsize'] = [10.,10]

mesh = [10,10,1]                                              # k-point mesh
a = 1.0                                                     # lattice constant    

def k_irr_BZ():    
    '''
    Calculates the k-vectors of the irreducable/full BZ:qdel 
    '''
    MAT = np.zeros((3,3))                                   # Matrix of reciprocal basis vectors               
    MAT[:,0] = np.array([1.,  1.0, 0.0])*2.*np.pi/a
    MAT[:,1] = np.array([1., -1.0, 0.0])*2.*np.pi/a
    MAT[:,2] = np.array([0.0, 0.0, 1.0])*2.*np.pi/a
        
    lattice = np.array([[1.0, 0.0, 0.0],                    # basis vectors of fcc lattice
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0]])*a  
                    
    positions = [[0.0, 0.0, 0.0]]                           # atomic basis in fractional coordinates             
       
    numbers= [1]                                      # atomic sorts (same digit for same atoms)   
        
    cell = (lattice, positions, numbers)
    
    print('spacegroup: ' +str(spglib.get_spacegroup(cell, symprec=1e-5)))
    print(spglib.get_symmetry(cell, symprec=1e-5))    
    
    # caclulatio of irr. BZ vectors + weights
    mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0,0,0])
    MAT_help = grid[np.unique(mapping)]/np.array(mesh, dtype=float)
    MAT_irr_BZ = np.zeros((np.size(MAT_help[:,0]),3))       
    for k in range(1,np.size(MAT_help[:,0])):
        MAT_irr_BZ[k,:] = MAT[:,0]*MAT_help[k,0] + MAT[:,1]*MAT_help[k,1] + MAT[:,2]*MAT_help[k,2] # transform from fractional to cartesian coordinates

    print("Number of kpoints: %d (irr BZ)" % len(np.unique(mapping)))
    num_kpoints = np.size(MAT_irr_BZ[:,0])

    weights = (np.unique(mapping,return_counts=True)[1])
    print("Number of kpoints: %d (full BZ, check of weights)" % weights.sum())           
           
    MAT_BZ_full = np.array(grid, dtype=float)
    for k in range(1,np.size(MAT_BZ_full[:,0])):
        MAT_BZ_full[k,:] = MAT[:,0]*MAT_BZ_full[k,0] + MAT[:,1]*MAT_BZ_full[k,1]+ MAT[:,2]*MAT_BZ_full[k,2]
    print("Number of kpoints: %d (full BZ)" % np.size(MAT_BZ_full[:,0]))
    
    file = open('k_BZ_irr.txt','w')
    for i in range(num_kpoints):
        for j in range(3):
            file.write("%s " % MAT_irr_BZ[i][j])
        file.write("\n")    
    file.close()
    
    file = open('k_weights_irr.txt','w')
    for i in range(num_kpoints):
        file.write("%s " % (weights[i]*1.0))
        file.write("\n")    
    file.close()

    file = open('k_weights_full.txt','w')
    for i in range(np.size(MAT_BZ_full[:,0])):
        file.write("%s " % 1.0)
        file.write("\n")    
    file.close()
    
    file = open('k_BZ_full.txt','w')
    for i in range(np.size(MAT_BZ_full[:,0])):
        for j in range(3):
            file.write("%s " % (MAT_BZ_full[i][j]/mesh[0]))
        file.write("\n")    
    file.close()
        
    return MAT_irr_BZ, MAT_BZ_full/(mesh[0]*1.0)



MAT_irr_BZ, MAT_BZ_full = k_irr_BZ()


############################################################################### PLOT DATA

import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import axes3d, Axes3D


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(MAT_irr_BZ[:,0], MAT_irr_BZ[:,1], MAT_irr_BZ[:,2], c="r", marker="x")
ax.scatter(MAT_BZ_full[:,0], MAT_BZ_full[:,1], MAT_BZ_full[:,2], c="k", marker=".")

ax.set_xlim(-6,+6)
ax.set_ylim(-6,+6)
plt.show()

