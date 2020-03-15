# Chern number calculation using the method of Fukui et al.

Author:
    [Damian Hofmann](damian.hofmann@mpsd.mpg.de)

Example code showing how to compute the Chern number (and, along the way, the Berry flux) for Hamiltonian on a discrete grid in 2d k-space.

This code uses the method described by [Fukui, Hatsugai, and Suzuki, J. Phys. Soc. Japan 74, 1674 (2005)](https://doi.org/10.1143/JPSJ.74.1674).

The code is written in [Julia](https://julialang.org) and has been tested with Julia version 0.6.3.

To run, obtain a working [Julia installation](https://julialang.org/downloads/) and execute

```shell
    $ julia fukui_berry_flux.jl
    Chern number for band #1: -1.0
    Chern number for band #2: 1.0
    Total Chern number: -0.0
```
