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