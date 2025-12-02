#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 2024

@author: felipe
"""

import numpy as np
import matplotlib.pyplot as plt
from functools import partial
from quadrature_lib import *

pn = lambda x, n: x**n
prim_pn = lambda x, n: (1/(n+1))*x**(n+1)

n_max = 6
a = 0.
b = 1.
dx = 0.1
N = int((b-a)/dx)

X=np.linspace(a, b, N)

methods = {'NC2': partial(newton_cotes, n=2),
           'NC3': partial(newton_cotes, n=3),
           'NC4': partial(newton_cotes, n=4),
           'NC5': partial(newton_cotes, n=5),
           'Ga1': partial(GLquadrature, n=1),
           'Ga2': partial(GLquadrature, n=2),
           'Ga3': partial(GLquadrature, n=3)}

Iex = np.array([prim_pn(b,k) - prim_pn(a,k) for k in range(n_max+1)])
Iapp = {}
error = {}

for m in methods.keys():
    quad_fun = methods[m]
    Iapp[m] = np.array([quad_fun(a,b,partial(pn,n=k)) for k in range(n_max+1)])
    error[m] = np.abs(Iapp[m] - Iex)

# un petite epsilon est nécessaire pour l'échelle logarithmique car log(0) n'est pas défini
eps = 1e-15
for m in methods.keys():
    plt.plot(error[m] + eps, '-o', label = m)
plt.grid()
plt.xlabel('ordre du polynôme')
plt.ylabel('erreur + eps')
plt.legend()
plt.yscale('log')
plt.savefig("degree_exactitude.png")
plt.savefig("degree_exactitude.pdf")
