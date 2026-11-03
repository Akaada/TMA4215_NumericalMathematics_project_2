import numpy as np

def RK4_step(f,x,v,h):
    
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
    B, DU = sys["B"], sys["DU"]

    def f(x, v):
        return np.cross(v, B(x)) - DU(x)

    N = int(round(T / h))
    X = np.empty((N + 1, 3))
    V = np.empty((N + 1, 3))
    X[0], V[0] = sys["x0"], sys["v0"]

    for n in range(N):
        X[n+1], V[n+1] = RK4_step(f, X[n], V[n], h)

    return X, V, h * np.arange(N + 1), "RK4"