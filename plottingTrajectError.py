import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def plotTraject(X,sys,methodName):
    nrows = 1
    if sys["name"] == "tokamak":
        nrows =2
    fig, axs  = plt.subplots(nrows=nrows,ncols=1,squeeze=False,figsize=(6, 5*nrows))


    axs[0].plot(X[:,0],X[:,1],lw = 0.5)
    axs[0].set_xlabel("$x_1$")
    axs[0].set_ylabel("$x_2$")
    axs[0].set_title(f"{sys['name']}: $(x_1, x_2)$ with {methodName} ")
    axs[0].set_aspect("equal")

    if sys["name"] == "tokamak":
        R = np.hypot(X[:,0],X[:,1])
        axs[1].plot(R,X[:,2],lw = 0.5)
        axs[1].set_xlabel("$R$"); axs[1].set_ylabel("$x_3$")
        axs[1].set_title("tokamak: $(R, x_3)$")

    fig.tight_layout()
    plt.show()


def errorPlot(X, times, Xref, sys,methodName ):
    err = np.linalg.norm(X - Xref, axis=1)
    plt.semilogy(times[1:], err[1:])   # skip n=0, because err = 0 and therefore log(0) = -inf
    plt.xlabel("$t_n$")
    plt.ylabel(r"$\|x^{(n)} - x_{\mathrm{ref}}(t_n)\|_2$")
    plt.title(f"{sys['name']}: error over time with {methodName}")
    plt.show()