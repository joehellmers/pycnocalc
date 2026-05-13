import csv
import os
import matplotlib.pyplot as plt

# --------------------------------------------------
# Files to compare
# --------------------------------------------------
files = [
    ("0.25", "researchdata/white_dwarf_benchmark_001_rstep_0p25.csv"),
    ("0.5",  "researchdata/white_dwarf_benchmark_001_rstep_0p5.csv"),
    ("1.0",  "researchdata/white_dwarf_benchmark_001_rstep_1p0.csv"),
]

# Benchmark point to pull from each CSV
target_rho = 5.0e9
target_temp = 1.0e8

# --------------------------------------------------
# Read one target row from a CSV
# --------------------------------------------------
def read_target_row(filename, rho_target, temp_target):
    with open(filename, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rho = float(row["rho"])
            temp = float(row["temp_k"])
            if abs(rho - rho_target) / rho_target < 1e-12 and abs(temp - temp_target) / temp_target < 1e-12:
                return {
                    "Rstep": float(row["Rstep"]),
                    "WKB_unscreened": float(row["WKB_unscreened"]),
                    "WKB_screened": float(row["WKB_screened"]),
                    "delta_WKB": float(row["delta_WKB"]),
                    "exp_delta_WKB": float(row["exp_delta_WKB"]),
                    "rate_ratio": float(row["rate_ratio"]),
                }
    return None

# --------------------------------------------------
# Load data
# --------------------------------------------------
data = []

for label, filename in files:
    if not os.path.exists(filename):
        print("Missing file:", filename)
        continue

    row = read_target_row(filename, target_rho, target_temp)
    if row is None:
        print("Target row not found in:", filename)
        continue

    data.append(row)

if len(data) < 2:
    raise RuntimeError("Need at least two data files to make a convergence plot.")

# Sort by Rstep ascending
data.sort(key=lambda x: x["Rstep"])

rsteps = [d["Rstep"] for d in data]
wkb_u = [d["WKB_unscreened"] for d in data]
wkb_s = [d["WKB_screened"] for d in data]
delta = [d["delta_WKB"] for d in data]
exp_delta = [d["exp_delta_WKB"] for d in data]
rate_ratio = [d["rate_ratio"] for d in data]

# --------------------------------------------------
# Relative error vs finest grid (smallest Rstep)
# --------------------------------------------------
ref_wkb_u = wkb_u[0]
ref_wkb_s = wkb_s[0]
ref_delta = delta[0]

err_wkb_u = [abs(x - ref_wkb_u) / abs(ref_wkb_u) * 100.0 for x in wkb_u]
err_wkb_s = [abs(x - ref_wkb_s) / abs(ref_wkb_s) * 100.0 for x in wkb_s]
err_delta = [abs(x - ref_delta) / abs(ref_delta) * 100.0 for x in delta]

# --------------------------------------------------
# Plot
# --------------------------------------------------
plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "legend.fontsize": 10,
})

fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)

fig.suptitle(
    r"Rstep Convergence for Mean-Field Screening" "\n"
    r"$^{12}$C + $^{12}$C, $\rho = 5\times10^9$ g/cm$^3$, $T = 10^8$ K",
    fontsize=16
)

# -----------------------------
# Panel 1: WKB values
# -----------------------------
ax = axes[0, 0]
ax.plot(rsteps, wkb_u, marker='o', linewidth=2, label='WKB unscreened')
ax.plot(rsteps, wkb_s, marker='o', linewidth=2, label='WKB screened')
ax.set_title("Raw WKB Values")
ax.set_xlabel("Rstep [fm]")
ax.set_ylabel("WKB")
ax.set_xticks(rsteps)
ax.grid(True, alpha=0.3)
ax.legend()

# -----------------------------
# Panel 2: delta WKB
# -----------------------------
ax = axes[0, 1]
ax.plot(rsteps, delta, marker='o', linewidth=2)
ax.set_title(r"Screening Effect: $\Delta$WKB")
ax.set_xlabel("Rstep [fm]")
ax.set_ylabel(r"$\Delta$WKB")
ax.set_xticks(rsteps)
ax.grid(True, alpha=0.3)

for x, y in zip(rsteps, delta):
    ax.annotate(f"{y:.4f}", (x, y), textcoords="offset points", xytext=(0, 6), ha='center', fontsize=9)

# -----------------------------
# Panel 3: enhancement
# -----------------------------
ax = axes[1, 0]
ax.plot(rsteps, exp_delta, marker='o', linewidth=2, label=r'$\exp(\Delta \mathrm{WKB})$')
ax.plot(rsteps, rate_ratio, marker='s', linewidth=2, linestyle='--', label='rate ratio')
ax.set_title("Screening Enhancement")
ax.set_xlabel("Rstep [fm]")
ax.set_ylabel("Enhancement")
ax.set_yscale("log")
ax.set_xticks(rsteps)
ax.grid(True, which='both', alpha=0.3)
ax.legend()

# -----------------------------
# Panel 4: convergence error
# -----------------------------
ax = axes[1, 1]
ax.plot(rsteps, err_wkb_u, marker='o', linewidth=2, label='WKB unscreened')
ax.plot(rsteps, err_wkb_s, marker='o', linewidth=2, label='WKB screened')
ax.plot(rsteps, err_delta, marker='o', linewidth=2, label=r'$\Delta$WKB')
ax.set_title("Relative Error vs Finest Step (0.25 fm)")
ax.set_xlabel("Rstep [fm]")
ax.set_ylabel("Percent error [%]")
ax.set_xticks(rsteps)
ax.grid(True, alpha=0.3)
ax.legend()

# Save
outfile = "researchdata/rstep_convergence_detailed.png"
plt.savefig(outfile, dpi=300, bbox_inches="tight")
plt.show()

print(f"Saved: {outfile}")
