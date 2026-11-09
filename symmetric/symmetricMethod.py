import sys as _sys   # aliased, since you use `sys` for the system dict
from pathlib import Path

# Parent folder of the current working directory:
_sys.path.insert(0, str(Path.cwd().parent))

# Or an explicit folder:
# _sys.path.insert(0, r"/full/path/to/folder")

from plottingTrajectError import findXref   # file name without .py

import numpy as np
from plottingTrajectError import findXref


def fpcd(y,j):
    return (y[j-2]-8*y[j-1]+8*y[j+1]-y[j+2])/12


def findF(j,x,Ax,sys,h):
    DA, DU = sys["DA"], sys["DU"]
    v = fpcd(x,j)/h
    Acd = fpcd(Ax,j)/h

    return DA(x[j]).T@v -Acd -DU(x[j]),v

def symmetricStep(i,x,Ax,F,h,sys,b0,b1,alpha):
    n = i-4
    change = (h**2)*(b1*F[n-1] +b0*F[n] + b1*F[n+1])
    return change -alpha[:8] @ x[n-4:n+4]

def symmetricMultistep(sys,h,T):
    a1 = -0.7
    a2 = 0.1
    a3 = 0.9

    b0 = (
        20*a1*a2*a3
        -4*(a1*a2+a1*a3+a2*a3)
        -28*(a1+a2+a3)
        -52)/3
    b1 = (
        2*a1*a2*a3
        +14*(a1*a2+a1*a3+a2*a3)
        +26*(a1 + a2 + a3)
        +38
    )/3
    alpha = np.array([1,-7/5,9/25,22/125,-34/125,22/125,9/25,-7/5,1])

    N = int(round(T / h))
    X  = np.empty((N + 3, 3))
    Ax = np.empty((N + 3, 3))
    F  = np.empty((N + 3, 3))
    V  = np.empty((N + 3, 3))
    Xs, Vs = findXref(sys, h*np.arange(8), rtol=1e-10, return_v=True)
    X[:8] = Xs
    V[:3] = Vs[:3]
    A = sys["A"]
    Ax[:8] = [A(xi) for xi in X[:8]]

    for j in (3, 4):
        F[j], V[j] = findF(j, X, Ax, sys, h)

    for i in range(8, N + 3):
        n = i - 4
        F[n+1], V[n+1] = findF(n+1, X, Ax, sys, h)
        X[i] = symmetricStep(i, X, Ax, F, h, sys, b0, b1, alpha)
        Ax[i] = A(X[i])

    V[N] = fpcd(X, N) / h

    return X[:N+1], V[:N+1],h * np.arange(N + 1), "Symmetric Multistep"



