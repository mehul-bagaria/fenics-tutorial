import numpy as np
import matplotlib.pyplot as plt

# Define function and its derivatives explicitly
def f(x):
    return np.sin(x)

def f1(x):
    return np.cos(x)

def f2(x):
    return -np.sin(x)

def f3(x):
    return -np.cos(x)

# Given parameters
x0 = np.pi / 4       # expansion point
h = np.pi / 8        # step
x_target = x0 + h    # point where we evaluate Taylor approximations

# Range for plotting: π/8 to π/2
x_vals = np.linspace(np.pi / 8, np.pi / 2, 400)
dx = x_vals - x0

# Taylor polynomials around x0 up to 3rd order
T1 = f(x0) + f1(x0) * dx
T2 = T1 + f2(x0) * dx**2 / 2.0
T3 = T2 + f3(x0) * dx**3 / 6.0

# True function and values at x_target
y_true = f(x_vals)
true_at_target = f(x_target)
T1_at_target = f(x0) + f1(x0) * (x_target - x0)
T2_at_target = T1_at_target + f2(x0) * (x_target - x0)**2 / 2.0
T3_at_target = T2_at_target + f3(x0) * (x_target - x0)**3 / 6.0

# Create the plot
plt.figure(figsize=(8, 6))

# Plot true function and Taylor approximations
plt.plot(x_vals, y_true,  linestyle="--", label="f(x) = sin(x)")
# plt.plot(x_vals, T1, linestyle="--", label="1st order Taylor")
# plt.plot(x_vals, T2, linestyle="--", label="2nd order Taylor")
# plt.plot(x_vals, T3, linestyle="--", label="3rd order Taylor")

# Mark x0 and x0 + h
plt.scatter([x0], [f(x0)], label="x0 = π/4")
plt.scatter([x_target], [true_at_target], label="true f(x0 + h)")

# Mark Taylor approximations at x0 + h
plt.scatter([x_target], [T1_at_target], label="T1(x0 + h)")
plt.scatter([x_target], [T2_at_target], label="T2(x0 + h)")
# plt.scatter([x_target], [T3_at_target], label="T3(x0 + h)")

plt.title("Taylor series approximations of sin(x) around x0 = π/4")
plt.xlabel("x")
plt.ylabel("f(x)")
plt.legend()
plt.grid(True)

plt.savefig("taylor_series_approximations.png")
plt.show()
