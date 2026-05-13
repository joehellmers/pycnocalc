import csv
import sys
import matplotlib.pyplot as plt


def read_data(filename):
    rho = []
    wkb_zeroT = []
    rate_zeroT = []
    turn1_zeroT = []
    turn2_zeroT = []
    turn1_zeroT_fm = []
    turn2_zeroT_fm = []

    with open(filename, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            rho.append(float(row["rho"]))
            wkb_zeroT.append(float(row["WKB_zeroT"]))
            rate_zeroT.append(float(row["rate_zeroT"]))
            turn1_zeroT.append(float(row["turn1_zeroT"]))
            turn2_zeroT.append(float(row["turn2_zeroT"]))
            turn1_zeroT_fm.append(float(row["turn1_zeroT_fm"]))
            turn2_zeroT_fm.append(float(row["turn2_zeroT_fm"]))

    return {
        "rho": rho,
        "wkb_zeroT": wkb_zeroT,
        "rate_zeroT": rate_zeroT,
        "turn1_zeroT": turn1_zeroT,
        "turn2_zeroT": turn2_zeroT,
        "turn1_zeroT_fm": turn1_zeroT_fm,
        "turn2_zeroT_fm": turn2_zeroT_fm,
    }


def make_plots(data, output_png):
    rho = data["rho"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    # WKB vs density
    ax = axes[0, 0]
    ax.plot(rho, data["wkb_zeroT"], marker="o")
    ax.set_xscale("log")
    ax.set_xlabel(r"Density $\rho$ [g/cm$^3$]")
    ax.set_ylabel("WKB")
    ax.set_title("Zero-T WKB vs Density")
    ax.grid(True)

    # Rate vs density
    ax = axes[0, 1]
    ax.plot(rho, data["rate_zeroT"], marker="o", color="darkred")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Density $\rho$ [g/cm$^3$]")
    ax.set_ylabel("Rate")
    ax.set_title("Zero-T Rate vs Density")
    ax.grid(True, which="both")

    # Turning point indices vs density
    ax = axes[1, 0]
    ax.plot(rho, data["turn1_zeroT"], marker="o", label="turn1 index")
    ax.plot(rho, data["turn2_zeroT"], marker="o", label="turn2 index")
    ax.set_xscale("log")
    ax.set_xlabel(r"Density $\rho$ [g/cm$^3$]")
    ax.set_ylabel("Turning point index")
    ax.set_title("Zero-T Turning Point Indices")
    ax.legend()
    ax.grid(True)

    # Turning point radii vs density
    ax = axes[1, 1]
    ax.plot(rho, data["turn1_zeroT_fm"], marker="o", label="turn1 [fm]")
    ax.plot(rho, data["turn2_zeroT_fm"], marker="o", label="turn2 [fm]")
    ax.set_xscale("log")
    ax.set_xlabel(r"Density $\rho$ [g/cm$^3$]")
    ax.set_ylabel("Turning point radius [fm]")
    ax.set_title("Zero-T Turning Point Radii")
    ax.legend()
    ax.grid(True)

    fig.tight_layout()
    fig.savefig(output_png, dpi=300)
    plt.show()


def main():
    if len(sys.argv) == 1:
        input_csv = "researchdata/turn_pt_zeroT_density_from_compare.csv"
        output_png = "researchdata/turn_pt_zeroT_density_from_compare.png"
    elif len(sys.argv) == 2:
        input_csv = sys.argv[1]
        output_png = "researchdata/turn_pt_zeroT_density_from_compare.png"
    elif len(sys.argv) == 3:
        input_csv = sys.argv[1]
        output_png = sys.argv[2]
    else:
        print("Usage: python3 plot_turn_pt_zeroT_density.py [input_csv] [output_png]")
        sys.exit(1)

    data = read_data(input_csv)
    make_plots(data, output_png)


if __name__ == "__main__":
    main()
