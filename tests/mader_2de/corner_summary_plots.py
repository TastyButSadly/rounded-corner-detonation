#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def grid_width(label: str) -> int:
    if "x" not in label:
        return 0
    head = label.split("x", 1)[0]
    return int(head) if head.isdigit() else 0


def read_rows(path: Path, label: str) -> list[dict[str, float | str]]:
    with path.open(newline="") as handle:
        rows = []
        for row in csv.DictReader(handle):
            rows.append({
                "grid": label,
                "R_corner_cm": float(row["R_corner_cm"]),
                "S_dead_cm2": float(row["S_dead_cm2"]),
                "S_hat": float(row["S_hat"]),
                "S_analysis_cm2": float(row["S_analysis_cm2"]),
                "S_box_cm2": float(row.get("S_box_cm2", "nan")),
                "arrived_fraction": float(row.get("arrived_fraction", "nan")),
            })
        return rows


def write_combined_csv(out_dir: Path, rows: list[dict[str, float | str]]) -> Path:
    path = out_dir / "summary_combined.csv"
    with path.open("w", newline="") as handle:
        fieldnames = ["grid", "R_corner_cm", "S_dead_cm2", "S_hat", "S_analysis_cm2", "S_box_cm2", "arrived_fraction"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    return path


def write_critical_csv(out_dir: Path, rows: list[dict[str, float | str]], epsilon: float) -> Path:
    path = out_dir / "critical_by_grid_summary.csv"
    labels = sorted({str(row["grid"]) for row in rows}, key=grid_width)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["grid", "epsilon", "R_cr_test_cm", "note"])
        writer.writeheader()
        for label in labels:
            data = sorted((row for row in rows if row["grid"] == label), key=lambda item: float(item["R_corner_cm"]))
            first = next((float(row["R_corner_cm"]) for row in data if float(row["S_hat"]) < epsilon), None)
            writer.writerow({
                "grid": label,
                "epsilon": epsilon,
                "R_cr_test_cm": f"{first:g}" if first is not None else f">{max(float(row['R_corner_cm']) for row in data):g}",
                "note": "first sampled radius below epsilon" if first is not None else "not found in sampled radii",
            })
    return path


def plot_radius_curves(out_dir: Path, rows: list[dict[str, float | str]], epsilon: float) -> None:
    labels = sorted({str(row["grid"]) for row in rows}, key=grid_width)
    for filename, y_key, y_label, title, min_width in (
        ("s_hat_by_radius.png", "S_hat", "S_dead / S_analysis [-]", "Relative dead-zone area vs radius", 0),
        ("s_hat_by_radius_fine.png", "S_hat", "S_dead / S_analysis [-]", "Relative dead-zone area vs radius, fine grids", 560),
        ("s_dead_by_radius_fine.png", "S_dead_cm2", "S_dead [cm^2]", "Dead-zone area vs radius, fine grids", 560),
    ):
        fig, axis = plt.subplots(figsize=(7.6, 4.8))
        for label in labels:
            if grid_width(label) < min_width:
                continue
            data = sorted((row for row in rows if row["grid"] == label), key=lambda item: float(item["R_corner_cm"]))
            axis.plot(
                [float(row["R_corner_cm"]) for row in data],
                [float(row[y_key]) for row in data],
                marker="o",
                linewidth=1.8,
                label=label,
            )
        if y_key == "S_hat":
            axis.axhline(epsilon, color="black", linestyle="--", linewidth=1.2, label=f"epsilon={epsilon:g}")
        axis.set_xlabel(r"$R_{\mathrm{corner}}$ [cm]")
        axis.set_ylabel(y_label)
        axis.set_title(title)
        axis.grid(True, alpha=0.3)
        axis.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(out_dir / filename, dpi=180)
        plt.close(fig)


def plot_arrived(out_dir: Path, rows: list[dict[str, float | str]]) -> None:
    labels = sorted({str(row["grid"]) for row in rows}, key=grid_width)
    fig, axis = plt.subplots(figsize=(7.6, 4.8))
    for label in labels:
        if grid_width(label) < 560:
            continue
        data = sorted((row for row in rows if row["grid"] == label), key=lambda item: float(item["R_corner_cm"]))
        axis.plot(
            [float(row["R_corner_cm"]) for row in data],
            [float(row["arrived_fraction"]) for row in data],
            marker="o",
            linewidth=1.8,
            label=label,
        )
    axis.set_xlabel("R_corner [cm]")
    axis.set_ylabel("S_analysis / S_box [-]")
    axis.set_title("Arrived fraction of downstream box, fine grids")
    axis.grid(True, alpha=0.3)
    axis.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "arrived_fraction_by_radius_fine.png", dpi=180)
    plt.close(fig)


def plot_grid_sensitivity(out_dir: Path, rows: list[dict[str, float | str]]) -> None:
    fig, axis = plt.subplots(figsize=(7.6, 4.8))
    for radius in sorted({float(row["R_corner_cm"]) for row in rows}):
        data = sorted((row for row in rows if abs(float(row["R_corner_cm"]) - radius) < 1e-12), key=lambda item: grid_width(str(item["grid"])))
        if len(data) < 2:
            continue
        axis.plot(
            [grid_width(str(row["grid"])) for row in data],
            [float(row["S_hat"]) for row in data],
            marker="o",
            linewidth=1.8,
            label=f"R={radius:g} cm",
        )
    axis.set_xlabel("Nx")
    axis.set_ylabel("S_dead / S_analysis [-]")
    axis.set_title("Grid sensitivity of dead-zone metric")
    axis.grid(True, alpha=0.3)
    axis.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "grid_sensitivity_summary.png", dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Create lightweight summary plots from rounded-corner study summaries.")
    parser.add_argument("--output-root", action="append", required=True)
    parser.add_argument("--label", action="append", required=True)
    parser.add_argument("--out-dir", default="tests/mader_2de/report/validation")
    parser.add_argument("--epsilon", type=float, default=0.01)
    args = parser.parse_args()

    if len(args.output_root) != len(args.label):
        raise SystemExit("--label must be supplied once per --output-root")

    rows: list[dict[str, float | str]] = []
    for root, label in zip(args.output_root, args.label):
        summary = Path(root).resolve() / "summary.csv"
        if not summary.exists():
            raise SystemExit(f"Missing summary: {summary}")
        rows.extend(read_rows(summary, label))

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"[write] {write_combined_csv(out_dir, rows)}")
    print(f"[write] {write_critical_csv(out_dir, rows, args.epsilon)}")
    plot_radius_curves(out_dir, rows, args.epsilon)
    plot_arrived(out_dir, rows)
    plot_grid_sensitivity(out_dir, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
