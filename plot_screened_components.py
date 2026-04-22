import csv
import matplotlib.pyplot as plt

r_vals = []
vc_vals = []
vs_vals = []
vn_vals = []

with open("researchdata/screened_components.csv", "r", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        r_vals.append(float(row["r_fm"]))
        vc_vals.append(float(row["Vcoulomb_MeV"]))
        vs_vals.append(float(row["Vscreen_MeV"]))
        vn_vals.append(float(row["Vnuclear_MeV"]))

plt.figure(figsize=(8, 5))
plt.plot(r_vals, vc_vals, label="V_Coulomb")
plt.plot(r_vals, vs_vals, label="V_screen")
plt.plot(r_vals, vn_vals, label="V_nuclear")

plt.xlabel("r [fm]")
plt.ylabel("Potential [MeV]")
plt.title("Coulomb, Screening, and Nuclear Potentials")
plt.xlim(0.5, 5.0)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
