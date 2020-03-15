#!/usr/bin/env python3

import numpy as np
from timeit import default_timer as timer

class Timer:
    """
    This class provides a timer that can be used in a
    with statement to time blocks of code.
    """
    def __init__(self, name, reference=None):
        self.name = name
        self.ref = reference
        self.last_delta = None

    def __enter__(self):
        self.start = timer()
        return self

    def __exit__(self, *args):
        delta = timer() - self.start
        self.last_delta = delta
        if self.ref:
            print(f"{self.name:30}: {delta:.4f} s (speedup: {self.ref/delta:.2f})")
        else:
            print(f"{self.name:30}: {delta:.4f} s")

from with_python import matmul_py
from with_ctypes import matmul_ctypes
from with_numba import matmul_numba1
from matmul_pybind import matmul_eigen, CustomMatrix

import pyximport; pyximport.install()
from with_cython import matmul_cy

def main():
    a = np.ones((200, 200))
    b = np.ones((200, 200))

    # This is slow, comment out if needed
    with Timer("Pure Python") as t:
        matmul_py(a, b)
    ref = t.last_delta

    with Timer("NumpPy", ref):
        np.dot(a, b)

    with Timer("ctypes", ref):
        matmul_ctypes(a, b)

    with Timer("cython", ref):
        matmul_cy(a, b)

    with Timer("numba", ref):
        matmul_numba1(a, b)

    a1 = CustomMatrix(200, 200)
    b1 = CustomMatrix(200, 200)
    with Timer("pybind11 / custom C++ class", ref):
        a1.dot(b1)

    with Timer("pybind11 / Eigen3", ref):
        matmul_eigen(a, b)

if __name__ == "__main__":
    main()
