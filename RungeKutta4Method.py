import numpy as np

def RK4_step(f,x,v,h):
    """
    Function to perform one step of RK4.
    Contained in its own function to make the main function less cluttered
    """

    # calculate the necesserary k-values
    kv1 = h*f(x,v)
    kx1 = h *v

    kv2 = h*f(x+0.5*kx1,v + 0.5*kv1)
    kx2 = h*(v+0.5*kv1)
    
    kv3 = h*f(x+0.5*kx2,v + 0.5*kv2)
    kx3 = h*(v+0.5*kv2)
    
    kv4 = h*f(x + kx3,v+kv3)
    kx4 = h*(v+kv3)



    v = v + (kv1 + 2*kv2 + 2*kv3 + kv4)/6
    x = x + (kx1 + 2*kx2 + 2*kx3 + kx4)/6
    return x,v


def RK4(sys, h, T):
    """The RK4 method that takes in the system as a dictionary, along with stepsize and end time"""
    # Get necessary functions related to the system
    B, DU = sys["B"], sys["DU"]

    # define the f function used by RK4
    def f(x, v):
        return np.cross(v, B(x)) - DU(x)


    N = int(round(T / h)) # Find number of steps
    X = np.empty((N + 1, 3)) # initalize X storage
    V = np.empty((N + 1, 3)) # Initalize V storage
    X[0], V[0] = sys["x0"], sys["v0"] # store initial values

    for n in range(N): # loop over number of steps
        X[n+1], V[n+1] = RK4_step(f, X[n], V[n], h) # do one step
    return X, V, h * np.arange(N + 1), "RK4" # Return the values necessary for future analysis of error and energy error. Along the name of the method to make life easier