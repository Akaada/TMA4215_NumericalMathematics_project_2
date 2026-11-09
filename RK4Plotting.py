from RungeKutta4Method import RK4
from plottingTrajectError import task5A
from Problem5B.problem5B import task5B


if __name__ == "__main__":  # plotter bare når filen kjøres direkte, ikke ved import
    # 5A: bane, feil over tid og energifeil (T = 1e5, h = 0.1)
    # save=True lagrer figurene som PNG i plots/
    task5A(RK4, "2d", save=True)
    task5A(RK4, "tokamak", save=True)

    # 5B: konvergensorden og work-precision (T = 1e3, h = 0.1, 0.05, 0.025)
    task5B([RK4], "2d", save=True)
    task5B([RK4], "tokamak", save=True)
