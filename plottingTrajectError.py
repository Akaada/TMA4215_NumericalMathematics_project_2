import numpy as np
import matplotlib.pyplot as plt

def plotTraject(X,sys):
    nrows = 1
    if sys["name"] == "tokamak":
        nrows =2
    fig, axs  = plt.subplots(nrows=nrows,ncols=1,squeeze=False)


    axs[0].plot(X[:,0],X[:,1])
    axs[0].set_xlabel("$x_1$")
    axs[0].set_ylabel("$x_2$")
    axs[0].set_title(f"{sys['name']}: $(x_1, x_2)$")
    axs[0].set_aspect("equal")

    if sys["name"] == "tokamak":
        R = np.hypot(X[:,0],X[:,1])
        axs[1].plot(R,X[:,2],lw = 0.5)
        axs[1].set_xlabel("$R$"); axs[1].set_ylabel("$x_3$")
        axs[1].set_title("tokamak: $(R, x_3)$")

    fig.tight_layout()
    plt.show()