import os
import numpy as np
from time import perf_counter
from create_systems import create_2d, create_tokamak
from plottingTrajectError import findXref
import matplotlib.pyplot as plt

# All 5B plots are saved in Problem5B/plots/
PLOT_DIR_5B = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")


def task5B(methods, system, T=1e3, h0=0.1, save=False):
    # Order of convergence and work-precision for a list of methods on one system ("2d" or "tokamak").
    # All methods are called as method(sys, h, T) and return (X, V, times, name),
    # so the function does not need to know which method it runs.
    # sys is a dict with the functions B, A, DA, U, DU, the initial values x0, v0 and the name ("2D"/"tokamak").
    # It is passed on to the methods, and sys['name'] is used in the printouts, titles and file names below
    sys = create_2d() if system == "2d" else create_tokamak()
    h = h0 / 2**np.arange(3)                 # h_i = h0 / 2^i: 0.1, 0.05, 0.025

    # Reference solution with solve_ivp (DOP853, rtol = atol = 1e-12). We only need the end point x_ref(T),
    # and T is the same for all runs, so it is computed once and reused for all methods and step sizes.
    XRef = findXref(sys,np.array([0.0,T]))
    print("XRef Done")

    # Left: error against step size (order). Right: runtime against error (work-precision).
    fig, axs = plt.subplots(1, 2, figsize=(12, 5))
    for method in methods:
        err = np.zeros(3)
        tau = np.zeros(3)
        for i, hi in enumerate(h):
            # The runtime is measured around the whole call, so the start-up (the half step backwards
            # for Boris, the starting values for the multistep method) is included
            start = perf_counter()
            X, V, times, name = method(sys, hi, T)
            tau[i] = perf_counter() - start
            # Global error at the end point: e_i = ||x_(N_i) - x_ref(T)||_2
            err[i] = np.linalg.norm(X[-1] - XRef[-1])
            print(f"{sys['name']:8s} {name:20s} h = {hi:.4f}: e = {err[i]:.3e}, tau = {tau[i]:.2f} s")

        # Observed order p_i = log2(e_{i-1} / e_i), since h is halved for each i:
        # e ≈ C h^p gives e_{i-1} / e_i ≈ 2^p
        orders = np.log2(err[:-1] / err[1:])
        print(f"{sys['name']:8s} {name:20s} observed order: " + ", ".join(f"{p:.2f}" for p in orders))

        # Log-log: e ≈ C h^p becomes a straight line with slope p
        axs[0].loglog(h, err, "o-", label=name)
        # Work-precision: tau ∝ e^(-1/p), so slope -1/p
        axs[1].loglog(err, tau, "o-", label=name)


    # Reference lines in black, so they are not confused with the methods.
    # They start at the first error of the last method in the list, and are only used to compare slopes.
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
        # File name with system and method names, e.g. 2D_RK4_boris_symmetricMultistep_order_workprecision.png
        names = "_".join(method.__name__ for method in methods)
        os.makedirs(PLOT_DIR_5B, exist_ok=True)
        fig.savefig(os.path.join(PLOT_DIR_5B, f"{sys['name']}_{names}_order_workprecision.png"),
                    dpi=150, bbox_inches="tight")
    plt.show()
    plt.close(fig)   # close the figure so the next plot is not drawn on top of it
