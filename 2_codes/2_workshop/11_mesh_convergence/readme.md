# Mesh Convergence Study for a First-Order ODE

## Problem Definition

We consider the ordinary differential equation:

$$
\frac{dy}{dx} = y, \qquad x \in [0,1],
$$

with the initial condition:

$$
y(0) = 1.
$$

The exact (analytical) solution is:

$$
y_{\text{exact}}(x) = e^{x}.
$$

This exact solution will be used to compute the error norms.

---

## Numerical Discretization

### Mesh Sizes

We perform the convergence study using four uniform mesh sizes:

$$
h \in \left\{\frac{1}{2},\ \frac{1}{4},\ \frac{1}{8},\ \frac{1}{16}\right\}.
$$

For each mesh size (h):

- the number of intervals is (N = 1/h),
- the grid points are (x_n = n h), for (n = 0,1,\dots,N),
- the numerical values will be denoted by (y_n).

### Forward Difference (Forward Euler) Approximation

We approximate the ODE using the explicit forward difference method:

$$
y_{n+1} = y_n + h \, y'(x_n, y_n),
$$

and for this problem, since (y' = y):

$$
y_{n+1} = y_n + h \, y_n = (1 + h) \, y_n.
$$

The initial value is:

$$
y_0 = 1.
$$

![](attachments/Pasted%20image%2020251126130414.png)

## Error Measures

We use the **L2 norm** to evaluate both:

1. the **local error** at a single point
2. the **global error** across the entire mesh

These error norms allow us to compute the observed convergence rate when the mesh is refined.



### Local L2 Error (Evaluated at x = 0.5)

For each mesh size (h), determine the grid node nearest to (x = 0.5):

- For (h = 1/2): node at 0.5 exists
- For (h = 1/4, 1/8, 1/16): node at 0.5 always exists because (0.5 = (0.5/h) \, h)

At that node:

$$
E_{\text{local}}(h) = \left| y_{\text{num}}(0.5) - e^{0.5} \right|.
$$

This is a pointwise (single-node) error but expressed as an L2 norm:

$$
|E_{\text{local}}(h)|_{2} = |E_{\text{local}}(h)|.
$$



### Global L2 Error (Evaluated Over All Nodes)

For each mesh, we compute:

$$
E_{\text{global}}(h) = \left( \sum_{n=0}^{N} \left( y_n - e^{x_n} \right)^2 \,  \right)^{1/2}
$$

This is the standard discrete L2 norm:

- error squared at every node,
- multiplied by the mesh size (h),
- square root taken at the end.

---
![](attachments/Pasted%20image%2020251126180410.png)
![](attachments/Pasted%20image%2020251126175503.png)
![](attachments/Pasted%20image%2020251126180359.png)

![](attachments/Pasted%20image%2020251126130517.png)