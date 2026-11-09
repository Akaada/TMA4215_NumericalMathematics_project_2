import os
import numpy as np
from time import perf_counter
from create_systems import create_2d, create_tokamak
from plottingTrajectError import findXref
import matplotlib.pyplot as plt

# Alle 5B-plott lagres i Problem5B/plots/
PLOT_DIR_5B = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")


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
            err[i] = np.linalg.norm(X[-1] - XRef[-1])
            print(f"{sys['name']:8s} {name:20s} h = {hi:.4f}: e = {err[i]:.3e}, tau = {tau[i]:.2f} s")

        # Observert orden p_i = log2(e_{i-1} / e_i), siden h halveres for hver i
        orders = np.log2(err[:-1] / err[1:])
        print(f"{sys['name']:8s} {name:20s} observed order: " + ", ".join(f"{p:.2f}" for p in orders))

        axs[0].loglog(h, err, "o-", label=name)
        axs[1].loglog(err, tau, "o-", label=name)


    # Referanselinjer i svart, så de ikke forveksles med metodene
    for p, style in ((2, "k--"), (4, "k:")):
        axs[0].loglog(h, err[0] * (h / h[0])**p, style, label=f"$h^{p}$")

    axs[0].set_xlabel("$h_i$"); axs[0].set_ylabel("$e_i$")
    axs[0].set_title(f"{sys['name']}: order verification")
    axs[1].set_xlabel("$e_i$"); axs[1].set_ylabel(r"$\tau_i$ [s]")
    axs[1].set_title(f"{sys['name']}: work–precision")
    for ax in axs:
        ax.legend()
    fig.tight_layout()
    if save:
        names = "_".join(method.__name__ for method in methods)
        os.makedirs(PLOT_DIR_5B, exist_ok=True)
        fig.savefig(os.path.join(PLOT_DIR_5B, f"{sys['name']}_{names}_order_workprecision.png"),
                    dpi=150, bbox_inches="tight")
    plt.show()
    plt.close(fig)


        







