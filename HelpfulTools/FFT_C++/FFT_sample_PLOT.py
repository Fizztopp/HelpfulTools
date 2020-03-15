#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 18:02:20 2017

@author: fizz
"""

## Fast Fourier Transform
import numpy as np
import matplotlib.pyplot as plt  

def func_test2(omega, eps, eta):
    return 1./(omega - eps + 1J*eta)

eps = 5.0
eta = 0.5

file = open('test_data.txt','r')
DATA = np.loadtxt(file)

fig = plt.figure(1, figsize=(10,10))

ax0 = fig.add_subplot(111)
OMEGA = np.linspace(-200,200,2000)
line1, = ax0.plot(OMEGA, func_test2(OMEGA, eps, eta).real, "b", label=r'real')
line2, = ax0.plot(OMEGA, func_test2(OMEGA, eps, eta).imag, "r", label=r'imag')
line3, = ax0.plot(OMEGA, DATA[:,0], 'b--', label=r'real (FFT)')
line4, = ax0.plot(OMEGA, DATA[:,1], "r--", label=r'imag (FFT)')
ax0.legend(fontsize=14)
ax0.set_xlim(-20,20)
plt.show()

