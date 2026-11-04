import numpy as np
from scipy.integrate import solve_ivp

def init_backward_step(sys, h):
    B, DU = sys["B"], sys["DU"]
    def f(t, y):
        x, v = y[:3], y[3:] 
        return np.concatenate([v,np.cross(v, B(x)) - DU(x)])
    y0 = np.concatenate([sys["x0"], sys["v0"]])
    sol = solve_ivp(f,(0.0, -h/2),y0,method = "DOP853",rtol=1e-10, atol=1e-12) #ja løsningen min bestod ikke antidoping (jk det er RK8)
    return sol.y[3:, -1]



    

def boris_step(x, v_half_back, h, sys):
    B, DU = sys["B"], sys["DU"]
    v_minus = v_half_back - 0.5*h*DU(x)          
    t = 0.5 * h * B(x)
    s = 2 * t / (1 + np.dot(t,t))
    v_prime = v_minus + np.cross(v_minus, t)
    v_plus = v_minus + np.cross(v_prime, s)    
    v_half_forward = v_plus - 0.5 * h * DU(x)      
    x_new = x + h * v_half_forward            
    return x_new, v_half_forward


def boris(sys, h, T):
    N = int(round(T / h))
    X = np.empty((N + 1, 3))
    V_half = np.empty((N + 2, 3))

    X[0] = sys["x0"]
    V_half[0] = init_backward_step(sys, h)

    for n in range(N):
        X[n + 1], V_half[n + 1] = boris_step(X[n], V_half[n], h, sys)

    _, V_half[N + 1] = boris_step(X[N], V_half[N], h, sys)

    V = 0.5 * (V_half[1:] + V_half[:-1])
    times = h * np.arange(N + 1)

    return X, V, times, "Boris"