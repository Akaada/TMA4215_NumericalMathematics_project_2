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
    """
    Calculating the five-point central difference centered in j
    Made so that we can give it y as an array
    """
    return (y[j-2]-8*y[j-1]+8*y[j+1]-y[j+2])/12


def findF(j,x,Ax,sys,h):
    """
    Finding F_j
    We assume that A(x) has been calculated already. We do this since A(x_j) is used 5 times each. So instead of calculating it 5 times we just store it
    """
    DA, DU = sys["DA"], sys["DU"] # Get the necessary functions from the system
    v = fpcd(x,j)/h # Calculate the v as the five-point central difference of x divided by step length h
    Acd = fpcd(Ax,j)/h # Calculate the A^(Delta) as the five-point central difference of A(x) divided by step length h

    return DA(x[j]).T@v -Acd -DU(x[j]),v # return the actual value for F_j

def symmetricStep(i,x,Ax,F,h,sys,b0,b1,alpha):
    """
    One step of the fourth-order symmetric multistep method
    We take in i as the last index of x we want to solve for
    """

    n = i-4 # We do an index shift so that we are now centered
    rhd = (h**2)*(b1*F[n-1] +b0*F[n] + b1*F[n+1]) # Calculate the right hand side first to make the return less cluttered
    return rhd -alpha[:8] @ x[n-4:n+4] # take the right hand side and subtract the first 8 parths of the lhd. The last alpha is 1, so no need for scaling

def symmetricMultistep(sys,h,T):
    """
    Solving a system with fourth order symmetric multistep
    Take in the system (sys) as a dictionary with necessary functions and initial values
    """

    # Initalizing the a-values
    a1 = -0.7
    a2 = 0.1
    a3 = 0.9


    # Calculating the b-values
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

    # Initalizing the alpha values
    alpha = np.array([1,-7/5,9/25,22/125,-34/125,22/125,9/25,-7/5,1])


    # Finding number of iterations we will need to do
    N = int(round(T / h))

    # Initializing matrixes for storing values
    # we go 2 beyond what we actually need since we use a five-point central difference to calculate v
    # And v is needed for energy calculations
    X  = np.empty((N + 3, 3))
    Ax = np.empty((N + 3, 3))
    F  = np.empty((N + 3, 3))
    V  = np.empty((N + 3, 3))

    # Use an exisisting function we have created that finds the first 7 additional values with a higher order RK method
    Xs, Vs = findXref(sys, h*np.arange(8), rtol=1e-10, return_v=True)

    # Store the new values
    X[:8] = Xs
    V[:3] = Vs[:3]

    A = sys["A"] # Take out the necessary function from the system dictionary
    Ax[:8] = [A(xi) for xi in X[:8]] # Calculating A(x) for the intial values

    for j in (3, 4):
        F[j], V[j] = findF(j, X, Ax, sys, h) # calculating the first necessary values for F and V

    for i in range(8, N + 3): # Loop over the number of iterations + 2
        n = i - 4 # Index shit necessary for calculating the new values for F and V
        F[n+1], V[n+1] = findF(n+1, X, Ax, sys, h) # Calculating F and V for the new step
        X[i] = symmetricStep(i, X, Ax, F, h, sys, b0, b1, alpha) # Do one step, and store the output
        Ax[i] = A(X[i]) # Store the evaluation of A in the new value of x

    V[N] = fpcd(X, N) / h # Calculate the last value of V


    # Return the postions, velocities, times, and name
    # Cut off since we went to steps too far
    return X[:N+1], V[:N+1],h * np.arange(N + 1), "Symmetric Multistep"



