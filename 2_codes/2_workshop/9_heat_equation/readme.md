# Heat Equation – README

## Introduction

The heat equation describes how temperature evolves inside a material over time. It is one of the most fundamental partial differential equations in physics, with applications in engineering, materials, meteorology, and numerical simulation.

In its simplest form, the heat equation expresses the balance between temporal change in temperature and spatial diffusion of heat.

## 1D Heat Equation

For a temperature field $T(x,t)$, the one-dimensional heat equation is:

$$
\frac{\partial T}{\partial t} = \alpha \frac{\partial^2 T}{\partial x^2}
$$

where

* $\alpha$ is the thermal diffusivity (material property),
* the left-hand side represents how temperature changes with time, and
* the right-hand side represents heat spreading due to conduction.

This equation governs how heat diffuses through a rod, wire, or any long body where heat transfer primarily occurs along one dimension.

## Intuition and Derivation (Simplified)

The equation arises from combining two ideas:

1. **Fourier's Law of Heat Conduction**
   Heat flux is proportional to the temperature gradient:
   
   $$
   q = -k \frac{\partial T}{\partial x}
   $$

2. **Energy Conservation in a Small Segment**
   Net heat entering a small element increases its internal energy:
   
   $$
   \rho c_p \frac{\partial T}{\partial t} = k \frac{\partial^2 T}{\partial x^2}
   $$

Dividing both sides by $\rho c_p$ gives the standard form with

$$
\alpha = \frac{k}{\rho c_p}
$$

## Example Solution Using Sine Functions
![](attachments/Pasted%20image%2020251125102420.png)

A simple and elegant solution of the heat equation is built using the sine function.
We use the fact that:

$$
\frac{d^2}{dx^2} \sin(x) = -\sin(x)
$$

Suppose temperature takes the form:

$$
T(x,t) = X(x) \cdot Y(t)
$$

Using separation of variables and applying the above identity:

* the spatial part gives:
  
  $$
  X''(x) = -\lambda X(x)
  $$
  
  where sine functions satisfy this naturally,

* the time part becomes:
  
  $$
  \frac{dY}{dt} = -\lambda \alpha \cdot Y(t)
  $$
  
  whose solution is an exponential decay:
  
  $$
  Y(t) = e^{-\lambda \alpha t}
  $$

Combining them:

$$
T(x,t) = \sin(x) \cdot e^{-\alpha t}
$$

This shows that temperature decays exponentially with time while maintaining its spatial sinusoidal shape.

![](attachments/Pasted%20image%2020251125102845.png)
## Another Simple Solution

If the temperature is linear in space (for example, $T(x) = ax + b$), its second derivative with respect to $x$ is zero.
Thus:

$$
\frac{\partial^2 T}{\partial x^2} = 0 \quad \Rightarrow \quad \frac{\partial T}{\partial t} = 0
$$

Meaning the temperature remains constant in time.
![](attachments/Pasted%20image%2020251125102911.png)
## Boundary Condition: No Heat Transfer

If the boundaries of the domain do not allow heat flow (insulated boundaries), then:

$$
-k \frac{\partial T}{\partial x} = 0 \quad \Rightarrow \quad \frac{\partial T}{\partial x} = 0
$$

This means the slope of the temperature curve is zero at the ends.
Such Neumann boundary conditions influence the form of valid solutions (e.g., cosine modes instead of sine modes).

![](attachments/Pasted%20image%2020251125104228.png)
## Reference Video

The following video provides a clear visual and mathematical explanation of the heat equation:
[https://www.youtube.com/watch?v=ToIXSwZ1pJU](https://www.youtube.com/watch?v=ToIXSwZ1pJU)
