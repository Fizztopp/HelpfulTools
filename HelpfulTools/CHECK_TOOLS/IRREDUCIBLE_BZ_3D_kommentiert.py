# Python programm to calculate the k-points of the full/irreducible Brilluoin zone. 
# Here for the example of an fcc lattice with a four atomic basis and 4 different atom sort
# (c) Gabriel E. Topp, gabriel.topp@mpsd.mpg.de

import numpy as np             
import spglib


mesh = [5,5,1]                                            # k-point mesh, 2d, third=1
a = 1.0                                                     # lattice constant    

def k_irr_BZ():    
    '''
    Calculates the k-vectors of the irreducable/full BZ:qdel  This was fcc, muss angepasst er
    '''
    MAT = np.zeros((3,3))                                   # Matrix of reciprocal basis vectors               
    MAT[:,0] = np.array([1., 1., 0.])*2.*np.pi/a
    MAT[:,1] = np.array([-1.,1., 0.])*2.*np.pi/a
    MAT[:,2] = np.array([ 0., 0.,1.])*2.*np.pi/a 			#Im prinziep irrelevant
        
    lattice = np.array([[1.0, 0.0, 0.0],                    # basis vectors of fcc lattice
                        [0.0, 0.1, 0.0],
                        [0.0, 0.0, 1.0]])*a  
    '''                
    positions = [[0.0, 0.0, 0.0],                           # atomic basis in fractional coordinates
                 [0.5, 0.0, 0.0],							# 2 Atome in der basis?
                 [0.0, 0.5, 0.0],
                 [0.0, 0.0, 0.5]]             
    '''
    
    positions = [[0.0, 0.0, 0.0]]
       
    #numbers= [1,2,3,4]                                      # atomic sorts (same digit for same atoms)   
    numbers=[1]
        
    cell = (lattice, positions, numbers)
    
    print("cell=", cell,'/n')
    '''
    cell= (array([[1. , 0. , 0. ],
				  [0. , 0.1, 0. ],
				  [0. , 0. , 1. ]]), [[0.0, 0.0, 0.0]], [1])
	'''
    
    print('spacegroup: ' +str(spglib.get_spacegroup(cell, symprec=1e-5))) #Gibt symmetriegruppen aus
    #print(spglib.get_symmetry(cell, symprec=1e-5))    
    
    # caclulation of irr. BZ vectors + weights  #Unter bestimmten bedingungen muss nicht die gesammte Reziproke Bz. 
    mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0,0,0])
    print("mapping=", mapping)
    print("grid=", grid)
    '''
    Irreducible k-points are obtained from a sampling mesh of k-points.
    mesh is given by three integers by array and specifies mesh numbers along reciprocal primitive axis. (mesh = [5,5,1])
    is_shift is given by the three integers by array. 
    When is_shift is set for each reciprocal primitive axis, the mesh is shifted along the axis in half of adjacent mesh points 
    irrespective of the mesh numbers. When the value is not 0, is_shift is set.
    
    grid gives the mesh points in fractional coordinates in reciprocal space.
    mapping gives mapping to the irreducible k-point indices that are obtained by np.unique(mapping)
    np.unique Finds the unique elements of an array.
    The grid point is accessed by grid[index].
    
    
    grid= [	[ 0  0  0]
			[ 1  0  0]
			[ 2  0  0]
			[-2  0  0]
			[-1  0  0]
			[ 0  1  0]
			[ 1  1  0]
			[ 2  1  0]
			[-2  1  0]
			[-1  1  0]
			[ 0  2  0]
			[ 1  2  0]
			[ 2  2  0]
			[-2  2  0]
			[-1  2  0]
			[ 0 -2  0]
			[ 1 -2  0]
			[ 2 -2  0]
			[-2 -2  0]
			[-1 -2  0]
			[ 0 -1  0]
			[ 1 -1  0]
			[ 2 -1  0]
			[-2 -1  0]
			[-1 -1  0]] das ist aber noch im real und nicht im Reziproken Raum
			
	mapping= [ 0  1  2  2  1  5  6  7  7  6 10 11 12 12 11 10 11 12 12 11  5  6  7  7 6] 
	sagt diese tabelle z.B, dass [-2  0  0] auf [2  0  0] abgebildet werden kann? Es sind 9 Zahlen und 9 kreuze in Rot
    '''
    
    MAT_help = grid[np.unique(mapping)]/np.array(mesh, dtype=float) #Teilen durch ein Array teilt x durch x usw. 
    
    
    print("np.unique(mapping)  ",np.unique(mapping))
    print("np.array(mesh, dtype=float)   ", np.array(mesh, dtype=float))
    print("grid[np.unique(mapping)]   ",grid[np.unique(mapping)])
    print("MAT_help", MAT_help)
    
    '''
    np.unique(mapping)   [ 0  1  2  5  6  7 10 11 12]
	np.array(mesh, dtype=float)    [5. 5. 1.]
	grid[np.unique(mapping)]    [	[0 0 0]
									[1 0 0]
									[2 0 0]
									[0 1 0]
									[1 1 0]
									[2 1 0]
									[0 2 0]
									[1 2 0]
									[2 2 0]]
								
	MAT_help [	[0.  0.  0. ]
				[0.2 0.  0. ]
				[0.4 0.  0. ]
				[0.  0.2 0. ]
				[0.2 0.2 0. ]
				[0.4 0.2 0. ]
				[0.  0.4 0. ]
				[0.2 0.4 0. ]
				[0.4 0.4 0. ]]  Warum geteilt wird verstehe ich nicht.

    '''
    
    MAT_irr_BZ = np.zeros((np.size(MAT_help[:,0]),3)) 
    print('np.size(MAT_help)   ',np.size(MAT_help))     			# np.size(MAT_help)    27
    print('np.size(MAT_help[:,0])   ',np.size(MAT_help[:,0]))     	# np.size(MAT_help[:,0])    9
    print('np.size(MAT_irr_BZ)   ',np.size(MAT_irr_BZ))				# np.size(MAT_irr_BZ)    27 
    print('np.shape(MAT_irr_BZ)   ',np.shape(MAT_irr_BZ))			# np.shape(MAT_irr_BZ)    (9, 3)

          
    for k in range(1,np.size(MAT_help[:,0])):
        MAT_irr_BZ[k,:] = MAT[:,0]*MAT_help[k,0] + MAT[:,1]*MAT_help[k,1] + MAT[:,2]*MAT_help[k,2] # transform from fractional to cartesian coordinates
        print('\n k=',k)
        print('MAT_irr_BZ[k,:]=',MAT_irr_BZ[k,:])
		# M  MAT_irr_BZ[k,:] ist eine Liste von Vektoren, welche den punkten der irr. Bz. im reziproken Raum entsprechen 		
		
    print("Number of kpoints: %d (irr BZ)" % len(np.unique(mapping)))
    num_kpoints_irr = np.size(MAT_irr_BZ[:,0])

    weights = (np.unique(mapping,return_counts=True)[1])
    print("Number of kpoints: %d (full BZ, check of weights)" % weights.sum())           
           
    MAT_BZ_full = np.array(grid, dtype=float)
    # Koordinaten in den Basiskoordinaten x mit xvektor, usw. 
    for k in range(1,np.size(MAT_BZ_full[:,0])):
        MAT_BZ_full[k,:] = MAT[:,0]*MAT_BZ_full[k,0] + MAT[:,1]*MAT_BZ_full[k,1]+ MAT[:,2]*MAT_BZ_full[k,2]	#Hier wird nicht geteilt
    print("Number of kpoints: %d (full BZ)" % np.size(MAT_BZ_full[:,0]))
    
    file = open('k_BZ_irr.txt','w')
    for i in range(num_kpoints_irr):
        for j in range(3):
            file.write("%s " % MAT_irr_BZ[i][j])
        file.write("\n")    
    file.close()
    
    file = open('k_weights_irr.txt','w')
    for i in range(num_kpoints_irr):
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
            file.write("%s " % (MAT_BZ_full[i][j]/np.array(mesh, dtype=float))) #Hier wird aber geteilt ???
        file.write("\n")    
    file.close()
        
    return MAT_irr_BZ, MAT_BZ_full/(mesh[0]*1.0)



MAT_irr_BZ, MAT_BZ_full = k_irr_BZ()


############################################################################### PLOT DATA

import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import axes3d, Axes3D


fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
ax = fig.add_subplot(111)
'''
ax.scatter(MAT_irr_BZ[:,0], MAT_irr_BZ[:,1], MAT_irr_BZ[:,2], c="r", marker="x")
ax.scatter(MAT_BZ_full[:,0], MAT_BZ_full[:,1], MAT_BZ_full[:,2], c="k", marker=".")
'''
ax.scatter(MAT_irr_BZ[:,0], MAT_irr_BZ[:,1], c="r", marker="x")
ax.scatter(MAT_BZ_full[:,0], MAT_BZ_full[:,1], c="k", marker=".")

plt.show()
