import csv
import sys
import matplotlib.pyplot as plt


def read_csv(filename):
    data = {
        "temp_k": [],
        "wkb_unscreened": [],
        "wkb_screened": [],
        "delta_wkb": [],
        "exp_delta_wkb": [],
        "turn1_unscreened": [],
        "turn2_unscreened": [],
        "turn1_screened": [],
        "turn2_screened": [],
        "turn1_unscreened_fm": [],
        "turn2_unscreened_fm": [],
        "turn1_screened_fm": [],
        "turn2_screened_fm": [],
    }

    with open(filename, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            data["temp_k"].append(float(row["temp_k"]))
            data["wkb_unscreened"].append(float(row["wkb_unscreened"]))
            data["wkb_screened"].append(float(row["wkb_screened"]))
            data["delta_wkb"].append(float(row["delta_wkb"]))
            data["exp_delta_wkb"].append(float(row["exp_delta_wkb"]))
            data["turn1_unscreened"].append(float(row["turn1_unscreened"]))
            data["turn2_unscreened"].append(float(row["turn2_unscreened"]))
            data["turn1_screened"].append(float(row["turn1_screened"]))
            data["turn2_screened"].append(float(row["turn2_screened"]))
            data["turn1_unscreened_fm"].append(float(row["turn1_unscreened_fm"]))
            data["turn2_unscreened_fm"].append(float(row["turn2_unscreened_fm"]))
            data["turn1_screened_fm"].append(float(row["turn1_screened_fm"]))
            data["turn2_screened_fm"].append(float(row["turn2_screened_fm"]))

    return data


def make_plots(data, output_png):
    x = data["temp_k"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    ax = axes[0, 0]
    ax.plot(x, data["wkb_unscreened"], marker="o", label="WKB unscreened")
    ax.plot(x, data["wkb_screened"], marker="o", label="WKB screened")
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("WKB")
    ax.set_title("WKB vs Temperature")
    ax.grid(True)
    ax.legend()

    ax = axes[0, 1]
    ax.plot(x, data["delta_wkb"], marker="o", color="purple")
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("delta WKB")
    ax.set_title("delta WKB vs Temperature")
    ax.grid(True)

    ax = axes[1, 0]
    ax.plot(x, data["exp_delta_wkb"], marker="o", color="darkred")
    ax.set_yscale("log")
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("exp(delta WKB)")
    ax.set_title("WKB Enhancement vs Temperature")
    ax.grid(True, which="both")

    ax = axes[1, 1]
    ax.plot(x, data["turn1_unscreened_fm"], marker="o", label="turn1 unscreened")
    ax.plot(x, data["turn2_unscreened_fm"], marker="o", label="turn2 unscreened")
    ax.plot(x, data["turn1_screened_fm"], marker="o", label="turn1 screened")
    ax.plot(x, data["turn2_screened_fm"], marker="o", label="turn2 screened")
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("Turning point radius [fm]")
    ax.set_title("Turning Points vs Temperature")
    ax.grid(True)
    ax.legend()

    fig.tight_layout()
    fig.savefig(output_png, dpi=300)
    plt.show()


def main():
    if len(sys.argv) == 1:
        input_csv = "researchdata/turn_pt_screened_compare.csv"
        output_png = "researchdata/turn_pt_screened_compare.png"
    elif len(sys.argv) == 2:
        input_csv = sys.argv[1]
        output_png = "researchdata/turn_pt_screened_compare.png"
    elif len(sys.argv) == 3:
        input_csv = sys.argv[1]
        output_png = sys.argv[2]
    else:
        print("Usage: python3 plot_turn_pt_screened_compare.py [input_csv] [output_png]")
        sys.exit(1)

    data = read_csv(input_csv)
    make_plots(data, output_png)


if __name__ == "__main__":
    main()
EOF
