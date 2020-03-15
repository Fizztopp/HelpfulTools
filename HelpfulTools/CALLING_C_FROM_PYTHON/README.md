# Calling C code from Python

This folder contains examples for several different ways of calling C or C++ routines from Python code.
(It has only been tested with Python 3.)

In order to run all examples, you should compile the C and C++ code by calling `make`.
You may need to install Python development headers on your system.
On Debian and Ubuntu, those should be contained in the package `python3-dev`.

In order to run the examples, you may need to install some Python modules, which you can do via `pip install cython numpy numba` (or by using the provided `requirements.txt` file via `pip install -r requirements.txt`).

In order to run all examples, call `./comparison.py`, which will output something like this:
```
$ pip install -r requirements.txt
$ make
$ ./comparison.py
Pure Python                   : 3.5430 s
NumpPy                        : 0.0004 s (speedup: 9297.89)
ctypes                        : 0.0398 s (speedup: 89.01)
cython                        : 0.0281 s (speedup: 126.07)
numba                         : 0.1776 s (speedup: 19.95)
pybind11 / custom C++ class   : 0.0090 s (speedup: 394.06)
pybind11 / Eigen3             : 0.0023 s (speedup: 1509.36)
```

Do not rely on these times (or the times you obtain from running this benachmark on your machine) as an accurate benchmark of the different approaches.
The examples here do not have the same level of optimization (actually, they are not really optimized at all, except probably for the functions calling the NumPy or Eigen3 code).

 Personally, I would recommend using [pybind11](https://github.com/pybind11/pybind11) for calling object-oriented C++ code from Python.
The C++ code one has to write to create the bindings is reasonably concise and the [documentation](https://pybind11.readthedocs.org) is quite good.
Furthermore, pybind11 provides automatic conversion of NumPy arrays to linear algebra classes from the [Eigen 3 framework](http://eigen.tuxfamily.org/index.php?title=Main_Page).
However, since it relies heavily on C++ template metaprogramming techniques, pybind11 can generate compiler error messages that are somewhat hard to read.
Also, it may not be the best option to bind C-style code that passes around lots of raw pointers.

Links for the other examples:

 * [Numba](http://numba.pydata.org/) does not call C/C++ code, but provides a just-in-time compiler to optimize numerical python code.
If the goal is to just speed-up a routine already implemented in Python, this might also provide a significant performance boost without much effort.

 * [Cython](https://cython.org) is an optimizing compiler that takes Python code (optinally with added data type annotations) and generates efficient C code.

 * [ctypes](https://docs.python.org/3.6/library/ctypes.html) is part of the Python standard library (so no packages need to be installed to use it), which allows calling C libraries from Python.
 At least in the example here, it seems to have some overhead, but if the goal is just to call into a C library from Python without high performance requirements, ctypes may be a quick option.

