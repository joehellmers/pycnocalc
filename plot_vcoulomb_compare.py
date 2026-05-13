import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------
# User inputs
# ---------------------------------------
Z1 = 6.0
Z2 = 6.0

# IMPORTANT:
# Use the same radii your Fortran code is using if you want an exact match.
# Replace these with your actual radius1 and radius2 values if you know them.
radius1 = 2.75   # fm
radius2 = 2.75   # fm

r_min = 0.1
r_max = 10.0
n_points = 1000

# ---------------------------------------
# Coulomb functions
# ---------------------------------------
def vcoulomb_finite(r, Z1, Z2, radius1, radius2):
    """
    Matches your Fortran Vcoulomb:
      if r < radius1 + radius2:
          V = ((3R^2 - r^2) * Z1*Z2*1.4397) / (2 R^3)
      else:
          V = (Z1*Z2*1.4397) / r
    """
    R = radius1 + radius2
    k = Z1 * Z2 * 1.4397

    if r < R:
        return ((3.0 * R**2 - r**2) * k) / (2.0 * R**3)
    else:
        return k / r


def vcoulomb_point(r, Z1, Z2):
    """
    Pure point-charge Coulomb potential: k/r
    """
    return (Z1 * Z2 * 1.4397) / r


# ---------------------------------------
# Generate data
# ---------------------------------------
r_values = np.linspace(r_min, r_max, n_points)
v_finite = np.array([vcoulomb_finite(r, Z1, Z2, radius1, radius2) for r in r_values])
v_point = np.array([vcoulomb_point(r, Z1, Z2) for r in r_values])

contact_radius = radius1 + radius2

# ---------------------------------------
# Plot
# ---------------------------------------
plt.figure(figsize=(8, 6))
plt.plot(r_values, v_finite, label='Finite-size Vcoulomb (code)')
plt.plot(r_values, v_point, linestyle='--', label='Point-charge Vcoulomb = Z1 Z2 1.4397 / r')
plt.axvline(contact_radius, linestyle=':', label=f'Contact radius = {contact_radius:.3f} fm')

plt.xlabel('r (fm)')
plt.ylabel('Potential (MeV)')
plt.title('Comparison of Finite-Size and Point-Charge Coulomb Potentials')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
