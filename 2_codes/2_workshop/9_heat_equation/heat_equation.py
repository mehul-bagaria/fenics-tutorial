"""
Simple 3D visualization of the 1D heat equation solution
T(x,t) = sin(x) * exp(-alpha * t)
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Parameters
alpha = 0.1  # Thermal diffusivity
x_max = 2 * np.pi  # Spatial domain [0, 2*pi]
t_max = 10.0  # Time domain [0, 10]

# Create mesh grid
x = np.linspace(0, x_max, 100)  # Space coordinate
t = np.linspace(0, t_max, 100)  # Time coordinate
X, T = np.meshgrid(x, t)

# Heat equation solution: T(x,t) = sin(x) * exp(-alpha * t)
Y = np.sin(X) * np.exp(-alpha * T)  # Temperature values

# Create 3D plot
fig = plt.figure(figsize=(12, 8))
ax = fig.add_subplot(111, projection='3d')

# Surface plot
surf = ax.plot_surface(X, Y, T, cmap='hot', alpha=0.8, edgecolor='none')

# Labels and title
ax.set_xlabel('X (Space Coordinate)', fontsize=12, labelpad=10)
ax.set_ylabel('Y (Temperature)', fontsize=12, labelpad=10)
ax.set_zlabel('Z (Time)', fontsize=12, labelpad=10)
ax.set_title('Heat Equation: T(x,t) = sin(x) · exp(-α·t)', fontsize=14, pad=20)

# Add colorbar
fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='Temperature')

# Adjust viewing angle
# elev: vertical rotation (positive tilts back)
# azim: horizontal rotation (0 = front view, looking along Y-axis)
ax.view_init(elev=20, azim=0)

plt.tight_layout()
plt.savefig("heat_equation_solution.png")
plt.show()