import numpy as np

def matmul_py(a, b):
    n1, m1 = a.shape
    n2, m2 = b.shape
    assert m1 == n2

    c = np.zeros((n1, m2))
    for i in range(n1):
        for j in range(m2):
            for k in range(m1):
                c[i,j] += a[i,k] * b[k,j]
    
    return c
