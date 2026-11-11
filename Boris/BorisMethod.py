import numpy as np
from scipy.integrate import solve_ivp

# The Boris method uses a staggered grid:
#   positions at the integer points   x_(n)     at t_n       = n*h
#   velocities at the half points     v_(n-1/2) at t_(n-1/2) = (n - 1/2)*h


def init_backward_step(sys, h):
    # The initial condition only gives x0 and v0 at t = 0, but the first Boris step needs v_(-1/2).
    # We find it by solving the first-order system y' = f(t, y), y = [x, v], backwards from t = 0 to t = -h/2.
    # sys is the dict from create_systems.py (create_2d or create_tokamak). We take out the functions
    # B(x)  and DU(x)  once, so f below can use them without looking them up every call
    B, DU = sys["B"], sys["DU"]
    def f(t, y):
        x, v = y[:3], y[3:]                                   # y has 6 components: 3 for x and 3 for v
        return np.concatenate([v,np.cross(v, B(x)) - DU(x)])  # [x', v'] = [v, v x B(x) - grad U(x)]
    y0 = np.concatenate([sys["x0"], sys["v0"]])   # initial values x0 and v0 from sys, stacked into one vector
    # DOP853 is an adaptive RK method of order 8. High accuracy here, since an error in the starting value
    # would otherwise propagate through the whole integration. It is only done once, so it is cheap.
    sol = solve_ivp(f,(0.0, -h/2),y0,method = "DOP853",rtol=1e-12, atol=1e-12) # yes, my solution failed the doping test (jk, it's RK8)
    return sol.y[3:, -1]   # last time point (t = -h/2), components 3-5 = v_(-1/2)


def boris_step(x, v_half_back, h, sys):
    # One Boris step: from x_(n) and v_(n-1/2) to x_(n+1) and v_(n+1/2).
    # The field functions B(x) and DU(x) from the system dict
    B, DU = sys["B"], sys["DU"]
    # Half a "kick" from the electric force -grad U
    v_minus = v_half_back - 0.5*h*DU(x)
    # Rotation of the velocity around B. t and s are vectors,
    t = 0.5 * h * B(x)
    s = 2 * t / (1 + np.dot(t,t))           # np.dot(t, t) = t^2
    v_prime = v_minus + np.cross(v_minus, t)
    v_plus = v_minus + np.cross(v_prime, s)
    # Second half kick from grad U gives the velocity at the next half point
    v_half_forward = v_plus - 0.5 * h * DU(x)
    # The position is updated with the velocity in the middle of the step
    x_new = x + h * v_half_forward
    return x_new, v_half_forward


def boris(sys, h, T):
    # The full integration from t = 0 to T with step size h.
    # Returns the same format as the other methods (X, V, times, name),
    # so the same plotting and error functions can be used for all of them.
    N = int(round(T / h))            # number of steps (round just in case)
    X = np.empty((N + 1, 3))         # x_(0), ..., x_(N)
    V_half = np.empty((N + 2, 3))    # V_half[k] = v_(k-1/2), i.e. v_(-1/2), ..., v_(N+1/2)

    X[0] = sys["x0"]                         # initial position from the system dict
    V_half[0] = init_backward_step(sys, h)

    for n in range(N):
        X[n + 1], V_half[n + 1] = boris_step(X[n], V_half[n], h, sys)

    # One extra step to get v_(N+1/2), which is needed for v_(N) below. The position is thrown away.
    _, V_half[N + 1] = boris_step(X[N], V_half[N], h, sys)

    # Velocities at the integer points (needed for the energy): v_(n) = (v_(n-1/2) + v_(n+1/2)) / 2
    V = 0.5 * (V_half[1:] + V_half[:-1])
    times = h * np.arange(N + 1)

    return X, V, times, "Boris"
