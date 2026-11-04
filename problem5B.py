import numpy as np
from time import perf_counter
from create_systems import create_2d, create_tokamak
from plottingTrajectError import findXref, showPlot
import matplotlib.pyplot as plt


def task5B(methods, system, T=1e3, h0=0.1, save=False):
    sys = create_2d() if system == "2d" else create_tokamak()
    h = h0 / 2**np.arange(3)                 
    XRef = findXref(sys,np.array([0.0,T]))
    print("XRef Done")

    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    for method in methods:
        err = np.zeros(3)
        tau = np.zeros(3)
        for i, hi in enumerate(h):
            start = perf_counter()
            X, V, times, name = method(sys, hi, T)
            tau[i] = perf_counter() - start
            print(f"i = {i} done")

            
            err[i] = np.linalg.norm(X[-1] - XRef[-1])   

        axs[0].loglog(h, err, "o-", label=name)
        axs[1].loglog(err, tau, "o-", label=name)



    for p in (2, 4):
        axs[0].loglog(h, err[0] * (h / h[0])**p, "--", label=f"$h^{p}$")

    axs[0].set_xlabel("$h_i$"); axs[0].set_ylabel("$e_i$")
    axs[0].set_title(f"{sys['name']}: order verification")
    axs[1].set_xlabel("$e_i$"); axs[1].set_ylabel(r"$\tau_i$ [s]")
    axs[1].set_title(f"{sys['name']}: work–precision")
    for ax in axs:
        ax.legend()
    fig.tight_layout()
    names = "_".join(method.__name__ for method in methods)
    showPlot(f"{sys['name']}_{names}_order_workprecision", save)


        







