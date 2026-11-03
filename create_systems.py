import numpy as np


def create_2d():
    def r(x):
        return np.hypot(x[0], x[1])

    def B(x):
        return np.array([0.0, 0.0, r(x)])

    def A(x):
        return r(x) / 3 * np.array([-x[1], x[0], 0.0])

    def DA(x):  
        R = r(x)
        return np.array([[-x[0]*x[1],    -(R**2 + x[1]**2), 0.0],
                         [R**2 + x[0]**2,  x[0]*x[1],       0.0],
                         [0.0,             0.0,             0.0]]) / (3*R)

    def U(x):
        return 1 / (100 * r(x))

    def DU(x):
        return -np.array([x[0], x[1], 0.0]) / (100 * r(x)**3)

    x0 = np.array([0.0, 1.0, 0.1])
    v0 = np.array([0.09, 0.05, 0.20])
    return dict(B=B, A=A, DA=DA, U=U, DU=DU, x0=x0, v0=v0, name="2D")


def create_tokamak():
    def U(x):
        return 0.0
    
    def DU(x):
        return np.zeros(3)


    def DA(x):
        x1, x2, x3 = x
        rho = np.hypot(x1, x2)
        return np.array([
            [0.0, 0.0, 0.0],
            [(rho - 1)/(2*rho), -x2/(2*rho*(rho + x1)), 0.0],
            [(x2*x3 - 2*x1)/(2*rho**2), -(x1*x3 + 2*x2)/(2*rho**2), -0.5*np.arctan(x2/x1)],
        ])

    def A(x):
        x1, x2, x3 = x
        rho = np.hypot(x1, x2)
        return np.array([
            0.0,
            0.5 * (x1 - np.log(rho + x1)),
            -0.5 * (x3 * np.arctan(x2 / x1) + np.log(rho**2)),
        ])

    def B(x):
        x1, x2, x3 = x
        rho = np.hypot(x1, x2)
        return np.array(
            [
                -(2*x2 + x1*x3)/(2*rho**2),
                (2*x1 - x2*x3)/(2*rho**2),
                (rho -1)/(2*rho)
            ]
        )

    x0 = np.array([1.2,0,0])
    v0 = np.array([0,4.816e-4,-2.059e-3])
    return dict(B=B, A=A, DA=DA, U=U, DU=DU, x0=x0, v0=v0, name="tokamak")