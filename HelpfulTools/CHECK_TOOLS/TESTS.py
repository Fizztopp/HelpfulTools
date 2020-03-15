#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 09:48:15 2020

@author: toppgabr
"""

import numpy as np  

A1 = np.array([[-1, 1+1j],[ 2-1j, 2]])
A2 = np.array([[-3, 2+2j],[ 3-3j, 4]])
print(np.dot(A1,A2))
print(np.dot(np.transpose(np.conj(A1)),A2))
print(np.dot(A1,np.transpose(np.conj(A2))))