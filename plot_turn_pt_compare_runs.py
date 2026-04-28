import csv
import sys
import matplotlib.pyplot as plt


def read_csv(filename):
    data = {
        "rho": [],
        "screen_temp": [],
        "Rstep": [],
        "partition": [],
        "nucIntType": [],
        "WKB_unscreened": [],
        "WKB_screened": [],
        "delta_WKB": [],
        "exp_delta_WKB": [],
        "turn1_unscreened": [],
        "turn2_unscreened": [],
        "turn1_screened": [],
        "turn2_screened": [],
    }

    with open(filename, "r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            for key in data:
                data[key].append(float(row[key]))

    return data


def plot_vs_temperature(data):
    x = data["screen_temp"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    ax = axes[0, 0]
    ax.plot(x, data["WKB_unscreened"], marker="o", label="WKB unscreened")
    ax.plot(x, data["WKB_screened"], marker="o", label="WKB screened")
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("WKB")
    ax.set_title("WKB vs Temperature")
    ax.grid(True)
    ax.legend()

    ax = axes[0, 1]
    ax.plot(x, data["delta_WKB"], marker="o", color="purple")
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("delta WKB")
    ax.set_title("delta WKB vs Temperature")
    ax.grid(True)

    ax = axes[1, 0]
    ax.plot(x, data["exp_delta_WKB"], marker="o", color="darkred")
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("exp(delta WKB)")
    ax.set_yscale("log")
    ax.set_title("WKB Enhancement vs Temperature")
    ax.grid(True, which="both")

    ax = axes[1, 1]
    ax.plot(x, data["turn1_unscreened"], marker="o", label="turn1 unscreened")
    ax.plot(x, data["turn2_unscreened"], marker="o", label="turn2 unscreened")
    ax.plot(x, data["turn1_screened"], marker="o", label="turn1 screened")
    ax.plot(x, data["turn2_screened"], marker="o", label="turn2 screened")
    ax.set_xlabel("Temperature [K]")
    ax.set_ylabel("Turning point index")
    ax.set_title("Turning Points vs Temperature")
    ax.grid(True)
    ax.legend()

    fig.tight_layout()
    plt.show()


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 plot_turn_pt_compare_runs.py <csv_file>")
        sys.exit(1)

    filename = sys.argv[1]
    data = read_csv(filename)
    plot_vs_temperature(data)


if __name__ == "__main__":
    main()
