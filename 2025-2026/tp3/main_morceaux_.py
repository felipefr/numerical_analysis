#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Nov 18 2024
@author: felipe
"""

from functools import partial
from quadrature_lib import *
import matplotlib.pyplot as plt
import numpy as np

c0 = 1.05
c1 = 1
c2 = 1
T1= 1
T2= 100

def f(t):
    return c0 + c1*np.exp(-t/T1) + c2*np.exp(-t/T2)

def prim_f(t):
    return # a compléter


def quad_morceaux(a,b,n,f,myquad):
    s = 0.;
    xi = np.linspace(a,b,n+1)
    for i in range(n):
        # à completer 
        
    return s


methods = {'NC2': partial(newton_cotes, n=2),
           'NC3': partial(newton_cotes, n=3),
           'NC4': partial(newton_cotes, n=4),
           'NC5': partial(newton_cotes, n=5),
           'Ga1': partial(GLquadrature, n=1),
           'Ga2': partial(GLquadrature, n=2),
           'Ga3': partial(GLquadrature, n=3)}


XY = np.loadtxt('Data0.txt');
X=XY[:,0]
Y=XY[:,1]
NX=len(X)
a, b = X.min() , X.max()

Yr = f(X)

plt.figure(1)
plt.plot(X,Y, label = 'mesures')
plt.plot(X,Yr, label = 'modele')
plt.legend()

# valeur exact
Iex = prim_f(b)-prim_f(a)


# à completer


