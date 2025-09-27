import numpy as np

def quad_NC2(a,b,f):
    s  = (b-a)*(f(a)+f(b))/2.
    return s

def quad_Ga1(a,b,f):
    x0 = (a+b)/2
    s  = (b-a)*f(x0)
    return s

# def quad_Ga2(a,b,f):
#     return s

# def quad_Ga3(a,b,f):
#     return s

# def quad_NC3(a,b,f):
#     return s

# def quad_NC4(a,b,f):
#     return s

# def quad_NC5(a,b,f):
#     return s

