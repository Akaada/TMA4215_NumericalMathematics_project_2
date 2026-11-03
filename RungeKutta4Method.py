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

