# -*- coding: utf-8 -*-
"""
Created on Fri Oct 12 13:34:16 2018

@author: toppgabr
"""

import matplotlib.pyplot as plt  
import numpy as np
import matplotlib.gridspec as gridspec
import matplotlib as mpl




file_ARPES= open('IPHOTO.txt','r')
IPHOTO = np.loadtxt(file_ARPES)
file_ARPES.close()


omega=np.linspace(-0.9,0.1,150) 
plt.plot(omega, IPHOTO)
plt.show()
    