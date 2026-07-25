from pathlib import Path
import sys

import matplotlib.pyplot as plt
import pandas as pd


# Read optional input/output paths supplied by rocket_sim. With no arguments,
# retain the original filenames so the script can still be run manually.
output_directory = Path(__file__).resolve().parent
if len(sys.argv) == 1:
    csv_path = output_directory / "flight_data.csv"
    png_path = output_directory / "flight_results.png"
elif len(sys.argv) == 3:
    csv_path = Path(sys.argv[1]).resolve()
    png_path = Path(sys.argv[2]).resolve()
else:
    raise SystemExit(
        "Usage: python plot_results.py [input_csv output_png]"
    )

data = pd.read_csv(csv_path)

# Create a two-by-two figure for the four requested flight plots.
figure, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
figure.suptitle("Rocket Flight Simulation Results", fontsize=16)

# Plot altitude throughout the flight.
axes[0, 0].plot(data["time"], data["altitude"], color="tab:blue")
axes[0, 0].set_title("Altitude vs Time")
axes[0, 0].set_xlabel("Time (s)")
axes[0, 0].set_ylabel("Altitude (m)")
axes[0, 0].grid(True, alpha=0.3)

# Plot vertical velocity throughout the flight.
axes[0, 1].plot(data["time"], data["velocity"], color="tab:orange")
axes[0, 1].set_title("Velocity vs Time")
axes[0, 1].set_xlabel("Time (s)")
axes[0, 1].set_ylabel("Velocity (m/s)")
axes[0, 1].grid(True, alpha=0.3)

# Overlay horizontal drift and fin correction using separate y-axes because
# they represent different units and may have different numerical scales.
drift_axis = axes[1, 0]
fin_axis = drift_axis.twinx()
drift_line = drift_axis.plot(
    data["time"], data["horizontal_drift"],
    color="tab:green", label="Horizontal drift"
)
fin_line = fin_axis.plot(
    data["time"], data["fin_correction_angle_deg"],
    color="tab:red", linestyle="--", label="Fin correction"
)
drift_axis.set_title("Horizontal Drift and Fin Correction")
drift_axis.set_xlabel("Time (s)")
drift_axis.set_ylabel("Horizontal Drift (m)", color="tab:green")
fin_axis.set_ylabel("Fin Correction (deg)", color="tab:red")
drift_axis.grid(True, alpha=0.3)
drift_axis.legend(
    drift_line + fin_line,
    [line.get_label() for line in drift_line + fin_line],
    loc="best",
)

# Plot the decreasing rocket mass during motor burn.
axes[1, 1].plot(data["time"], data["mass"], color="tab:purple")
axes[1, 1].set_title("Mass vs Time")
axes[1, 1].set_xlabel("Time (s)")
axes[1, 1].set_ylabel("Mass (kg)")
axes[1, 1].grid(True, alpha=0.3)

# Save a presentation-quality image beside the input CSV and close the figure.
figure.savefig(png_path, dpi=300)
plt.close(figure)
