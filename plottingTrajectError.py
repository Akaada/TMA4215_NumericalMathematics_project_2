import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from create_systems import create_2d, create_tokamak


 ### 
 #dette er bare for å lagre plot underveis, kan fjernes når vi lager notebook, da bytter vi bare showplot func med plt.show
# Figurer lagres i mappen plots/ i rotmappen til prosjektet
PLOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "plots")


def showPlot(filename, save=False):
    """Viser gjeldende figur. Hvis save=True lagres den også som plots/<filename>.png."""
    if save:
        os.makedirs(PLOT_DIR, exist_ok=True)
        plt.savefig(os.path.join(PLOT_DIR, f"{filename}.png"), dpi=150, bbox_inches="tight")
    plt.show()
    plt.close()   # lukk figuren så neste plott ikke tegnes oppå denne

###
def plotTraject(X, sys, methodName, save=False):
    nrows = 2 if sys["name"] == "tokamak" else 1
    fig, axs = plt.subplots(nrows=nrows, ncols=1, squeeze=False, figsize=(6, 5*nrows))
    axs = axs[:, 0]   # squeeze=False gives a 2D array; take the column of axes

    axs[0].plot(X[:, 0], X[:, 1], lw=0.5)
    axs[0].set_xlabel("$x_1$")
    axs[0].set_ylabel("$x_2$")
    axs[0].set_title(f"{sys['name']}: $(x_1, x_2)$ with {methodName}")
    axs[0].set_aspect("equal")

    if sys["name"] == "tokamak":
        R = np.hypot(X[:, 0], X[:, 1])
        axs[1].plot(R, X[:, 2], lw=0.5)
        axs[1].set_xlabel("$R$")
        axs[1].set_ylabel("$x_3$")
        axs[1].set_title(f"tokamak: $(R, x_3)$ with {methodName}")

    fig.tight_layout()
    showPlot(f"{sys['name']}_{methodName}_trajectory", save)


def errorPlot(X, times, Xref, sys, methodName, save=False):
    err = np.linalg.norm(X - Xref, axis=1)
    plt.figure()   # ny figur, ellers tegnes plottet i forrige figur
    plt.semilogy(times[1:], err[1:])   # skip n=0: err = 0 there, and log(0) = -inf
    plt.xlabel("$t_n$")
    plt.ylabel(r"$\|x^{(n)} - x_{\mathrm{ref}}(t_n)\|_2$")
    plt.title(f"{sys['name']}: error over time with {methodName}")
    showPlot(f"{sys['name']}_{methodName}_error", save)


def findXref(sys, times,rtol = 1e-12):
    B, DU = sys["B"], sys["DU"]

    def f(t, y):
        x, v = y[:3], y[3:]
        return np.concatenate([v, np.cross(v, B(x)) - DU(x)])

    y0 = np.concatenate([sys["x0"], sys["v0"]])
    solution = solve_ivp(f, (times[0], times[-1]), y0, method="DOP853",
                         t_eval=times, rtol=rtol, atol=1e-12)
    if not solution.success:
        print("Reference solution failed:", solution.message)

    # solution.y has shape (6, N+1); rows 0-2 are positions. Transpose to (N+1, 3).
    return solution.y[:3].T


def E(X, V, U):
    kinetic = 0.5 * np.sum(V**2, axis=1)
    potential = np.array([U(x) for x in X])
    return potential + kinetic


def plotEnergyError(X, V, U, times, sys, methodName, save=False):
    energy = E(X, V, U)
    error = np.abs(energy - energy[0])
    plt.figure()   # ny figur, ellers tegnes plottet i forrige figur
    plt.semilogy(times[1:], error[1:])
    plt.xlabel("$t_n$")
    plt.ylabel("$e_n$")
    plt.title(f"{sys['name']}: energy error with {methodName}")
    showPlot(f"{sys['name']}_{methodName}_energy", save)


def task5A(method, system, Xref=None, T=1e5, h=0.1, save=False):
    sys = create_2d() if system == "2d" else create_tokamak()

    X, V, times, methodName = method(sys, h, T)
    if Xref is None:
        Xref = findXref(sys, times)

    plotTraject(X, sys, methodName, save)
    errorPlot(X, times, Xref, sys, methodName, save)
    plotEnergyError(X, V, sys["U"], times, sys, methodName, save)


