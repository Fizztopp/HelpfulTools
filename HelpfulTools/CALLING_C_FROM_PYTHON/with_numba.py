import numpy as np
from timeit import default_timer as timer
import numba

class Timer:
    def __init__(self, name):
        self.name = name

    def __enter__(self):
        self.start = timer()
        return self

    def __exit__(self, *args):
        delta = timer() - self.start
        print(f"{self.name}: {delta:.4f} seconds")
        

@numba.jit
def matmul_numba1(a, b):
    n1, m1 = a.shape
    n2, m2 = b.shape
    assert m1 == n2

    c = np.zeros((n1, m2))
    for i in range(n1):
        for j in range(m2):
            for k in range(m1):
                c[i,j] += a[i,k] * b[k,j]
    
    return c


def main():
    a = np.random.rand(300, 5000)
    b = np.random.rand(5000, 300)

    with Timer("Numba"):
        c1 = matmul_numba1(a, b)
        
    with Timer("NumPy"):
        c2 = np.dot(a, b)

    assert np.allclose(c1, c2)

if __name__ == "__main__":
    main()
