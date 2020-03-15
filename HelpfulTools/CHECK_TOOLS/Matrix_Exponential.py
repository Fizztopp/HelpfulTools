# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 16:56:34 2019

@author: toppgabr
"""

import numpy as np  
from scipy.linalg import expm, norm               
import spglib
import matplotlib as mpl
import matplotlib.pyplot as plt 
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import axes3d, Axes3D

dt = 0.1

MATRIX = np.random.rand(2,2)+1J*np.random.rand(2,2)
MATRIX = 0.5*(MATRIX+np.transpose(np.conj(MATRIX)))
#print(np.linalg.eigh(MATRIX))

MAT_DIAG = np.array([[np.linalg.eigh(MATRIX)[0][0],0],[0,np.linalg.eigh(MATRIX)[0][1]]])
MAT_EXP = np.array([[np.exp(-1J*dt*np.linalg.eigh(MATRIX)[0][0]),0],[0,np.exp(-1J*dt*np.linalg.eigh(MATRIX)[0][1])]])
S = np.linalg.eigh(MATRIX)[1]
SCT  = np.transpose(np.conj(S))


print(MATRIX)
print("1")
print(np.dot(S,np.dot(MAT_DIAG,SCT)))
print("2")
print(expm(-1J*MATRIX*dt))
print("3")
print(np.dot(S,np.dot(MAT_EXP,SCT)))
print("4")
print(np.dot(S,np.dot(MAT_EXP,SCT))-expm(-1J*MATRIX*dt))
print("5")
print(1j*np.log(np.linalg.eig(expm(-1J*MATRIX*dt))[0]))                        # NO Hermitian matrix!!! --> do not use eigh()!!!
print("6")
print(np.log(np.linalg.eigh(expm(MATRIX*dt))[0]))                              # NO Hermitian matrix!!! --> do not use eigh()!!!
print("7")
print(np.linalg.eig(expm(-1J*MATRIX*dt))[1])
print("8")
print(np.linalg.eig(1J*MATRIX)[1])
print("")


S = np.linalg.eigh(np.array([[0,-2+1j],[-2-1j,0]]))[1]
SCT  = np.transpose(np.conj(S))
print("")
print("")
print(np.linalg.eigh(np.array([[0,-2+1j],[-2-1j,0]]))[0])
print("")
print(np.dot(S,np.dot(np.array([[0,-2+1j],[-2-1j,0]]),SCT)))                   # wrong! 
print(np.dot(SCT,np.dot(np.array([[0,-2+1j],[-2-1j,0]]),S)))                   # eigenvectors are the columns of S, H_d = S.H.S^(-1)!!!

#Chiral Symmetrie s3.H.s3 = -H for mass=0 (Graphene)#
mass = 0.1
kx = 0.1
ky = 0.2

s0 = np.eye(2)
s1 = np.array([[0., 1.],[1., 0.]])
s2 = np.array([[0.,-1j],[1j, 0.]])
s3 = np.array([[1., 0.],[0.,-1.]])

HMAT = kx*s1 + ky*s2 + mass*s3

print(np.dot(np.conj(S[:,0]),np.dot(np.array([[0,-2+1j],[-2-1j,0]]),S[:,0])))