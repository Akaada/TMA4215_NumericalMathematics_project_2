import os
import sys as pysys   # pysys, siden vi bruker `sys` for system-dicten

# Rotmappen og metodemappene må ligge på søkestien for at importene skal virke
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
for folder in (ROOT, os.path.join(ROOT, "Boris"), os.path.join(ROOT, "symmetric")):
    pysys.path.insert(0, folder)

from RungeKutta4Method import RK4
from BorisMethod import boris
from symmetricMethod import symmetricMultistep
from problem5B import task5B   


if __name__ == "__main__":
    # 5B for alle tre metodene i samme plott, lagres i Problem5B/plots/
    methods = [RK4, boris, symmetricMultistep]
    task5B(methods, "2d", save=True)
    task5B(methods, "tokamak", save=True)
