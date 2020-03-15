import numpy as np
cimport numpy as np

def matmul_cy(np.ndarray[np.double_t, ndim=2] a, np.ndarray[np.double_t, ndim=2] b):
    cdef int n1 = a.shape[0]
    cdef int m1 = a.shape[1]
    cdef int n2 = b.shape[0]
    cdef int m2 = b.shape[1]
    assert m1 == n2

    cdef np.ndarray[np.double_t, ndim=2] c
    c = np.zeros((n1, m2))
    
    cdef int i, j, k
    for i in range(n1):
        for j in range(m2):
            for k in range(m1):
                c[i,j] += a[i,k] * b[k,j]
    
    return c
