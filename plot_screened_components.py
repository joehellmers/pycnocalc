import csv
import matplotlib.pyplot as plt

r_vals = []
vc_vals = []
vs_vals = []
vc_screened_vals = []
vn_vals = []
vtot_unscreened_vals = []
vtot_screened_vals = []

with open("researchdata/screened_components.csv", "r", newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        r_vals.append(float(row["r_fm"]))
        vc_vals.append(float(row["Vcoulomb_MeV"]))
        vs_vals.append(float(row["Vscreen_MeV"]))
        vc_screened_vals.append(float(row["Vcoulomb_screened_MeV"]))
        vn_vals.append(float(row["Vnuclear_MeV"]))
        vtot_unscreened_vals.append(float(row["Vtotal_unscreened_MeV"]))
        vtot_screened_vals.append(float(row["Vtotal_screened_MeV"]))

plt.figure(figsize=(8, 5))

plt.plot(r_vals, vc_vals, label="V_Coulomb")
plt.plot(r_vals, vn_vals, label="V_nuclear")
plt.plot(r_vals, vtot_unscreened_vals, label="V_total unscreened", linewidth=2)
plt.plot(r_vals, vtot_screened_vals, label="V_total screened", linewidth=2)

# Optional: keep these if you want to show the screening piece explicitly
# plt.plot(r_vals, vs_vals, label="V_screen", linestyle="--")
# plt.plot(r_vals, vc_screened_vals, label="V_Coulomb screened", linestyle=":")

plt.xlabel("r [fm]")
plt.ylabel("Potential [MeV]")
plt.title("Unscreened vs Screened Total Potential")
plt.xlim(0.5, 5.0)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
