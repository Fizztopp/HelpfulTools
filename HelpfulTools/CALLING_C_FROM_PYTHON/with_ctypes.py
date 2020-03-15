import ctypes
import numpy as np

lib = ctypes.cdll.LoadLibrary('./libmatmul.so')

def matmul_ctypes(a, b):
    n1, m1 = a.shape
    n2, m2 = b.shape
    assert m1 == n2

    c = np.ndarray((n1, m2))

    ap = a.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    bp = b.ctypes.data_as(ctypes.POINTER(ctypes.c_double))
    cp = c.ctypes.data_as(ctypes.POINTER(ctypes.c_double))

    lib.matmul(ap, bp, cp, n1, m1, m2)

    return c
