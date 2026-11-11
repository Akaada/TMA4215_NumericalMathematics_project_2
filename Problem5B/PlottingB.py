import os
import sys as pysys   # pysys, since we use `sys` for the system dict

# The root folder and the method folders must be on the search path for the imports to work
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for folder in (ROOT, os.path.join(ROOT, "Boris"), os.path.join(ROOT, "symmetric")):
    pysys.path.insert(0, folder)

from RungeKutta4Method import RK4
from BorisMethod import boris
from symmetricMethod import symmetricMultistep
from problem5B import task5B


if __name__ == "__main__":
    # 5B for all three methods in the same plot, saved in Problem5B/plots/
    # The methods are passed in as functions, task5B calls them with (sys, h, T)
    methods = [RK4, boris, symmetricMultistep]
    task5B(methods, "2d", save=True)
    task5B(methods, "tokamak", save=True)
