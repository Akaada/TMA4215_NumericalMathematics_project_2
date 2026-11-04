import os
import sys as pysys   # kalt pysys for å ikke forveksles med system-dicten sys i funksjonene

# Fellesfilene (create_systems, plottingTrajectError, problem5B) ligger i mappen over Boris
# dette gjør så denne filen kan se de
pysys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from BorisMethod import boris
from plottingTrajectError import task5A
from problem5B import task5B


if __name__ == "__main__":
    # 5A: bane, feil over tid og energifeil (T = 1e5, h = 0.1)
    # save=True lagrer figurene som PNG i plots/ (i rotmappen)
    task5A(boris, "2d", save=True)
    task5A(boris, "tokamak", save=True)

    # 5B: konvergensorden og work-precision (T = 1e3, h = 0.1, 0.05, 0.025)
    task5B([boris], "2d", save=True)
    task5B([boris], "tokamak", save=True)
