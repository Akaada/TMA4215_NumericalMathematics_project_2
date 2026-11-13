import numpy as np


def create_2d():
    """Function to create the necessary functions and intial values for the 2D system. Stored in a dictionary"""


    def r(x): # a small function used to calculate the hypotenuse. Not used in later code
        return np.hypot(x[0], x[1])

    def B(x): # The B function
        return np.array([0.0, 0.0, r(x)])

    def A(x): # The A function
        return r(x) / 3 * np.array([-x[1], x[0], 0.0])

    def DA(x):  # The jacobian of the A function. Calculated by hand. Derivation is included in the notebook
        R = r(x)
        return np.array([[-x[0]*x[1],    -(R**2 + x[1]**2), 0.0],
                         [R**2 + x[0]**2,  x[0]*x[1],       0.0],
                         [0.0,             0.0,             0.0]]) / (3*R)

    def U(x): # The U function
        return 1 / (100 * r(x))

    def DU(x): # The derviative of the U function. Given in the task
        return -np.array([x[0], x[1], 0.0]) / (100 * r(x)**3)


    # inital values
    x0 = np.array([0.0, 1.0, 0.1])
    v0 = np.array([0.09, 0.05, 0.20])
    return dict(B=B, A=A, DA=DA, U=U, DU=DU, x0=x0, v0=v0, name="2D") # A dictionary to make life easier, along with name of system for ease of plotting


def create_tokamak():
    """Function to create the necessary functions and intial values for the tokamak system. Stored in a dictionary"""

    def U(x): # The U function. Defaulted to zero. Created since we will use functions that need a function for U. 
        # Is a possibility to have them stored as numbers, and then instead use if statements that check if U is a function or numbers. But we found this easier
        return 0.0
    
    def DU(x):
        # See U(x)
        return np.zeros(3) # return a matrix since the functions expect a matrix.


    def DA(x): # The jacobian of the A function. Calculated by hand. Derivation is included in the notebook
        x1, x2, x3 = x
        rho = np.hypot(x1, x2)
        return np.array([
            [0.0, 0.0, 0.0],
            [(rho - 1)/(2*rho), -x2/(2*rho*(rho + x1)), 0.0],
            [(x2*x3 - 2*x1)/(2*rho**2), -(x1*x3 + 2*x2)/(2*rho**2), -0.5*np.arctan2(x2, x1)],
        ])

    def A(x): # The A function
        x1, x2, x3 = x
        rho = np.hypot(x1, x2)
        return np.array([
            0.0,
            0.5 * (x1 - np.log(rho + x1)),
            -0.5 * (x3 * np.arctan2(x2, x1) + np.log(rho**2)),
        ])

    def B(x): # The B function
        x1, x2, x3 = x
        rho = np.hypot(x1, x2)
        return np.array(
            [
                -(2*x2 + x1*x3)/(2*rho**2),
                (2*x1 - x2*x3)/(2*rho**2),
                (rho -1)/(2*rho)
            ]
        )

    # inital values
    x0 = np.array([1.2,0,0])
    v0 = np.array([0,4.816e-4,-2.059e-3])
    return dict(B=B, A=A, DA=DA, U=U, DU=DU, x0=x0, v0=v0, name="tokamak") # A dictionary to make life easier, along with name of system for ease of plotting