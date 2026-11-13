## Structure

1. Introduction
    1. Motivation
    2. Equations of motion
    3. Systems
    
3. Numerical methods
    1. Boris method
    2. Classical Runge–Kutta 4 (RK4)
    3. Fourth-order symmetric multistep method
4. Implementation
5. Long-time simulations
    1. Trajectories
    2. Global error
    3. Energy conservation
6. Order verification and efficiency
    1. Order of convergence
    2. Work–precision diagrams
7. Exploratory experiment




# Introduction 

## Motivation 
Simulating charged particles in magnetic fields requires integrating velocity and position over very long time periods. Standard methods like the Runge-Kutta method are accurate over short times, but do not preserve the geometric structure of the problem. Therefore small errors add up over long periods, systematically. This causes the energy to slowly diverge from its true value and the simulation becomes unphysical.  
In this project we compare RK4 with two methods designed for long-time simulations of charged particles: the Boris method and a fourth-order symmetric multistep method. 
All the methods and problems discussed are as in the given project description.

## Problems 
We use the notation $R := \sqrt{x_1^2 + x_2^2}$.

### 2.1 2D charged particle 
The magnetic field, vector potential and electrostatic potential are

$$
B(x) = \begin{pmatrix} 0 \\ 0 \\ R \end{pmatrix}, \qquad
A(x) = \frac{1}{3}\begin{pmatrix} -x_2 R \\ x_1 R \\ 0 \end{pmatrix},
$$

$$
U(x) = \frac{1}{100\,R}, \qquad
\nabla U(x) = -\frac{1}{100\,R^3}\begin{pmatrix} x_1 \\ x_2 \\ 0 \end{pmatrix}.
$$

The initial values are

$$
x(0) = \begin{pmatrix} 0.0 \\ 1.0 \\ 0.1 \end{pmatrix}, \qquad
v(0) = \begin{pmatrix} 0.09 \\ 0.05 \\ 0.20 \end{pmatrix}.
$$

### 2.2 Tokamak magnetic field 
The magnetic field and vector potential are

$$
B(x) = \begin{pmatrix}
-\dfrac{2x_2 + x_1 x_3}{2R^2} \\[2ex]
\dfrac{2x_1 - x_2 x_3}{2R^2} \\[2ex]
\dfrac{R - 1}{2R}
\end{pmatrix}, \qquad
A(x) = \begin{pmatrix}
0 \\[1ex]
\dfrac{1}{2}\Big(x_1 - \ln(R + x_1)\Big) \\[2ex]
-\dfrac{1}{2}\Big(x_3 \arctan\!\big(\tfrac{x_2}{x_1}\big) + \ln(x_1^2 + x_2^2)\Big)
\end{pmatrix},
$$

and there is no electrostatic potential, $U(x) = \nabla U(x) = 0$.
The energy (3) then reduces to $E = \frac{1}{2}\|v\|_2^2$. The initial
values are

$$
x(0) = \begin{pmatrix} 1.2 \\ 0 \\ 0 \end{pmatrix}, \qquad
v(0) = \begin{pmatrix} 0 \\ 4.816 \cdot 10^{-4} \\ -2.059 \cdot 10^{-3} \end{pmatrix}.
$$

Note that $A(x)$ is singular for $x_2 = 0$, $x_1 \le 0$ due to the
logarithm. With the given initial values, the particle does not cross
this region.


## Systems
**2D charged particle.** The magnetic field points along the $x_3$-axis
and increases in strength with the distance $R$ from the axis, while the
electrostatic force $-\nabla U$ points radially outwards. The particle
gyrates around the field lines while slowly drifting around the
$x_3$-axis. Since neither $B$ nor $\nabla U$ has a force component in
the $x_3$-direction, $v_3$ stays constant, and all interesting dynamics
takes place in the $(x_1, x_2)$-plane.

**Tokamak magnetic field.** The field consists of a toroidal part,
circling the $x_3$-axis with magnetic strength proportional to $1/R$, and a
poloidal part, circling the magnetic axis $R = 1$, $x_3 = 0$. Together
they produce a helical field, with lines winding around a torus. The particle
starts on the outer side of the magnetic axis with small parallel
velocity, and is reflected back and fourth by the stronger field on the inner side, tracing out a so-called banana
orbit in the $(R, x_3)$-plane.

# Numerical methods


### Boris method
The Boris method is a second-order method on a staggered grid.
The position is computed at integer time steps $t_n$, while the velocity
is computed at half steps $t_{n+1/2}$. Given $x^{(n)}$, $v^{(n-1/2)}$ and
step size $h$, one step consists of a half kick from the electric field,
a rotation of the velocity around the magnetic field, a second half kick
and finally an update of the position:

$$
\begin{aligned}
v^- &= v^{(n-1/2)} - \frac{h}{2}\nabla U(x^{(n)}), \qquad
t = \frac{h}{2} B(x^{(n)}), \\
v' &= v^- + v^- \times t, \\
v^+ &= v^- + v' \times s, \qquad s = \frac{2t}{1 + t \cdot t}, \\
v^{(n+1/2)} &= v^+ - \frac{h}{2}\nabla U(x^{(n)}), \\
x^{(n+1)} &= x^{(n)} + h\, v^{(n+1/2)}.
\end{aligned}
$$

Since only $x^{(0)}$ and $v^{(0)}$ are given, we initialise
$v^{(-1/2)}$ by integrating (1) half a step backwards in time with
`solve_ivp` (DOP853,
`rtol=1e-12`, `atol=1e-12`). The velocity at integer time steps is found by averaging,

$$
v^{(n)} = \frac{1}{2}\left(v^{(n+1/2)} + v^{(n-1/2)}\right).
$$

### Classical Runge–Kutta 4 (RK4)
RK4 is applied to the first-order system (1). With
$f(x, v) = v \times B(x) - \nabla U(x)$, one step from
$(x^{(n)}, v^{(n)})$ with step size $h$ is given by

$$
\begin{aligned}
k_1^v &= h\, f\big(x^{(n)},\, v^{(n)}\big),
& k_1^x &= h\, v^{(n)}, \\
k_2^v &= h\, f\big(x^{(n)} + \tfrac{1}{2}k_1^x,\, v^{(n)} + \tfrac{1}{2}k_1^v\big),
& k_2^x &= h\big(v^{(n)} + \tfrac{1}{2}k_1^v\big), \\
k_3^v &= h\, f\big(x^{(n)} + \tfrac{1}{2}k_2^x,\, v^{(n)} + \tfrac{1}{2}k_2^v\big),
& k_3^x &= h\big(v^{(n)} + \tfrac{1}{2}k_2^v\big), \\
k_4^v &= h\, f\big(x^{(n)} + k_3^x,\, v^{(n)} + k_3^v\big),
& k_4^x &= h\big(v^{(n)} + k_3^v\big),
\end{aligned}
$$

$$
\begin{aligned}
v^{(n+1)} &= v^{(n)} + \frac{1}{6}\left(k_1^v + 2k_2^v + 2k_3^v + k_4^v\right), \\
x^{(n+1)} &= x^{(n)} + \frac{1}{6}\left(k_1^x + 2k_2^x + 2k_3^x + k_4^x\right).
\end{aligned}
$$

### Fourth-order symmetric multistep method 
This is an explicit 9-point method for the second-order formulation (2).
Instead of the magnetic field $B(x)$, it uses the vector potential
$A(x)$ and its Jacobian $A'(x)$, which we compute analytically.
Let the five-point central difference centred at $y^{(j)}$ be

$$
\Delta y^{(j)} := \frac{1}{12}\left(y^{(j-2)} - 8y^{(j-1)} + 8y^{(j+1)} - y^{(j+2)}\right).
$$

With $x^{(j)} \approx x(t_j)$, $t_j = t_0 + jh$, we define

$$
v^{(j)} := \frac{\Delta x^{(j)}}{h}, \qquad
A_j^\Delta := \frac{\Delta A(x^{(j)})}{h}, \qquad
F_j := A'(x^{(j)})^\top v^{(j)} - A_j^\Delta - \nabla U(x^{(j)}).
$$

The method is then given by

$$
\sum_{j=-4}^{4} \alpha_j\, x^{(n+j)} = h^2\left(\beta_1 F_{n-1} + \beta_0 F_n + \beta_1 F_{n+1}\right)
$$

with coefficients

$$
(\alpha_{-4}, \dots, \alpha_4) = \left(1, -\tfrac{7}{5}, \tfrac{9}{25}, \tfrac{22}{125}, -\tfrac{34}{125}, \tfrac{22}{125}, \tfrac{9}{25}, -\tfrac{7}{5}, 1\right),
$$

$$
\begin{aligned}
\beta_0 &= \tfrac{1}{3}\big(20a_1a_2a_3 - 4(a_1a_2 + a_1a_3 + a_2a_3) - 28(a_1 + a_2 + a_3) - 52\big), \\
\beta_1 &= \tfrac{1}{3}\big(2a_1a_2a_3 + 14(a_1a_2 + a_1a_3 + a_2a_3) + 26(a_1 + a_2 + a_3) + 38\big),
\end{aligned}
$$

where $a_1 = -0.7$, $a_2 = 0.1$, $a_3 = 0.9$ as in [2]. This gives
$\beta_0 = -19.74$ and $\beta_1 = 12.378$.

Since $\alpha_4 = 1$ and $F_{n+1}$ only depends on
$x^{(n-1)}, \dots, x^{(n+3)}$, the method is explicit, and the new
position is found as

$$
x^{(n+4)} = h^2\left(\beta_1 F_{n-1} + \beta_0 F_n + \beta_1 F_{n+1}\right) - \sum_{j=-4}^{3} \alpha_j\, x^{(n+j)}.
$$

The method needs the eight starting values $x^{(0)}, \dots, x^{(7)}$,
which we generate by solving with `solve_ivp` (DOP853,
`rtol=1e-12`, `atol=1e-12`).


### Jacobian of the 2D vector potential

We use the same convention $(A'(x))_{ij} = \dfrac{\partial A_i}{\partial x_j}$ as in tokamak. With $r = \sqrt{x_1^2 + x_2^2}$, the vector potential is

$$
A(x) = \frac{r}{3}\begin{pmatrix} -x_2 \\ x_1 \\ 0 \end{pmatrix}, \qquad
A_1 = -\tfrac13 x_2 r, \qquad A_2 = \tfrac13 x_1 r, \qquad A_3 = 0.
$$



$$
\frac{\partial r}{\partial x_1} = \frac{x_1}{r}, \qquad
\frac{\partial r}{\partial x_2} = \frac{x_2}{r}, \qquad
\frac{\partial r}{\partial x_3} = 0.
$$

No component depends on $x_3$, therefore the third column is zero.  Since $A_3 = 0$, the third row is also zero. The remaining entries follow from the product rule.

**Row 1:**

$$
\frac{\partial A_1}{\partial x_1} = -\frac13 x_2 \frac{x_1}{r} = -\frac{x_1 x_2}{3r}, \qquad
\frac{\partial A_1}{\partial x_2} = -\frac13\Big(r + x_2\frac{x_2}{r}\Big) = -\frac{r^2 + x_2^2}{3r}.
$$

**Row 2:**

$$
\frac{\partial A_2}{\partial x_1} = \frac13\Big(r + x_1\frac{x_1}{r}\Big) = \frac{r^2 + x_1^2}{3r}, \qquad
\frac{\partial A_2}{\partial x_2} = \frac13 x_1 \frac{x_2}{r} = \frac{x_1 x_2}{3r}.
$$

**Result**

$$
A'(x) = \frac{1}{3r}
\begin{pmatrix}
-x_1 x_2 & -(r^2 + x_2^2) & 0 \\
r^2 + x_1^2 & x_1 x_2 & 0 \\
0 & 0 & 0
\end{pmatrix}
$$

The 3D matrix explains a 2D system, i.e only the upper left corner has values. This makes sense considering we are in a 2D system. 

### Jacobian of the Tokamak vector potential

We use the convention $(A'(x))_{ij} = \dfrac{\partial A_i}{\partial x_j}$, so row $i$ is the gradient of the component $A_i$. We define $\rho = \sqrt{x_1^2 + x_2^2}$ to make notation easier. The vector potential in the Tokamak system is defined as

$$
A_1 = 0, \qquad
A_2 = \tfrac12\big(x_1 - \ln(\rho + x_1)\big), \qquad
A_3 = -\tfrac12\Big(x_3 \arctan\tfrac{x_2}{x_1} + \ln(\rho^2)\Big).
$$

**Helping derivatives**
We calculate some derivatives that will be used troughout the derivation. Making it easier to do the actual calculation.
$$
\frac{\partial \rho}{\partial x_1} = \frac{x_1}{\rho}, \qquad
\frac{\partial \rho}{\partial x_2} = \frac{x_2}{\rho}, \qquad
\frac{\partial}{\partial x_1}\arctan\frac{x_2}{x_1} = -\frac{x_2}{\rho^2}, \qquad
\frac{\partial}{\partial x_2}\arctan\frac{x_2}{x_1} = \frac{x_1}{\rho^2},
$$

$$
\frac{\partial}{\partial x_1}\ln(\rho^2) = \frac{2x_1}{\rho^2}, \qquad
\frac{\partial}{\partial x_2}\ln(\rho^2) = \frac{2x_2}{\rho^2}.
$$

**Row 1:** $A_1 = 0$, making all entries on row 1 be zero.

**Row 2:** $A_2$ does not depend on $x_3$, so $\partial A_2/\partial x_3 = 0$. For $x_1$:

$$
\frac{\partial}{\partial x_1}\ln(\rho + x_1) = \frac{x_1/\rho + 1}{\rho + x_1} = \frac{(x_1 + \rho)/\rho}{\rho + x_1} = \frac{1}{\rho}
\quad\Rightarrow\quad
\frac{\partial A_2}{\partial x_1} = \frac12\Big(1 - \frac{1}{\rho}\Big) = \frac{\rho - 1}{2\rho}.
$$

For $x_2$:

$$
\frac{\partial A_2}{\partial x_2} = -\frac12 \cdot \frac{x_2/\rho}{\rho + x_1} = -\frac{x_2}{2\rho(\rho + x_1)}.
$$

**Row 3:**

$$
\frac{\partial A_3}{\partial x_1} = -\frac12\Big(-\frac{x_2 x_3}{\rho^2} + \frac{2x_1}{\rho^2}\Big) = \frac{x_2 x_3 - 2x_1}{2\rho^2},
$$

$$
\frac{\partial A_3}{\partial x_2} = -\frac12\Big(\frac{x_1 x_3}{\rho^2} + \frac{2x_2}{\rho^2}\Big) = -\frac{x_1 x_3 + 2x_2}{2\rho^2},
$$

$$
\frac{\partial A_3}{\partial x_3} = -\frac12 \arctan\frac{x_2}{x_1}.
$$

**Result**

$$
A'(x) =
\begin{pmatrix}
0 & 0 & 0 \\[4pt]
\dfrac{\rho - 1}{2\rho} & -\dfrac{x_2}{2\rho(\rho + x_1)} & 0 \\[10pt]
\dfrac{x_2 x_3 - 2x_1}{2\rho^2} & -\dfrac{x_1 x_3 + 2x_2}{2\rho^2} & -\dfrac12 \arctan\dfrac{x_2}{x_1}
\end{pmatrix}
$$

