# args is only used when f has additional arguments
# for the basic case, f = f(y,t) only.
# For the vectorial case, y0 should be a list. 
# If y0 is a list, y is matrix. If y0 is scalar, y is a vector 

def euler_exp(f, y0, t, args=()):
    n = len(t)
    y = np.zeros((n, len(y0))) if type(y0) is list else np.zeros(n)
    y[0] = y0
    for i in range(n - 1):
        y[i+1] = y[i] + (t[i+1] - t[i]) * f(y[i], t[i], *args)
    
    return y

def rungekutta4(f, y0, t, args=()):
    n = len(t)
    y = np.zeros((n, len(y0))) if type(y0) is list else np.zeros(n)
    y[0] = y0
    for i in range(n - 1):
        h = t[i+1] - t[i]
        k1 = f(y[i], t[i], *args)
        k2 = f(y[i] + k1 * h / 2., t[i] + h / 2., *args)
        k3 = f(y[i] + k2 * h / 2., t[i] + h / 2., *args)
        k4 = f(y[i] + k3 * h, t[i] + h, *args)
        y[i+1] = y[i] + (h / 6.) * (k1 + 2*k2 + 2*k3 + k4)
    return y

# Verlet's method: also knows leap-frog (saute-mouton en français)
def verlet(f, y0, t, args = ()):
    n = len(t)
    y = np.zeros((n, 2)) # only valid for second order
    y[0,:] = y0 # correspond to time t0
    
    for i in range(n-1):
        h = t[i+1] - t[i]
        aold = f(y[i], t[i], *args)[1]
        y[i+1,0] = y[i,0] + h*y[i,1] + 0.5*h**2*aold # f_1 is the acceleration
        anew = f(y[i+1], t[i+1], *args)[1]
        y[i+1,1] = y[i,1] + 0.5*h*(anew + aold)
    
    return y

# F for the pendulum
def pend(y, t, k1, k2):
    return np.array([y[1], -k1*y[1] - k2*np.sin(y[0])])

def energy_kinetics(y, m,l):
    return 0.5*m*(l*y[1])**2

def energy_potential(y, m, g, l): 
    return m*g*l*(1-np.cos(y[0]))
