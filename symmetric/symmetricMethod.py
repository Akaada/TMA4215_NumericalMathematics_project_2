import numpy as np
from scipy.integrate import solve_ivp
from plottingTrajectError import findXref


def fpcd(y,j):
    return (y[j-2]-8*y[j-1]+8*y[j+1]-y[j+2])/12

def startVals(sys,i=7):
    return findXref(sys,np.array([0,i]),rtol=1e-10)

def findF(j,x,Ax,DA,DU,h):
    v = fpcd(x,j)/h
    Acd = fpcd(Ax,j)/h

    return DA(x[j]).T@v -Acd -DU(x[j])

def symmetricStep(x,i):
    j = i-4
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
    )
    alpha = np.array([1,-7/5,9/25,22/125,-34/125,22/125,9/25,-7/5,1])

    
