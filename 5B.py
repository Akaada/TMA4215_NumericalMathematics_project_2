import numpy as np
from time import perf_counter
from create_systems import create_2d, create_tokamak
from plottingTrajectError import findXref


def findError(X,Xref):
    return  np.linalg.norm(X[-1] - Xref[-1])



def task5B(method,system,XRef = None,T = 1e3):
    h0 = 0.1
    sys = create_2d() if system == "2d" else create_tokamak()


    for i in [0,1,2]:
        h = h0/(2**i)
        

        if XRef == None:
            XRefi = findXref(sys)


        







