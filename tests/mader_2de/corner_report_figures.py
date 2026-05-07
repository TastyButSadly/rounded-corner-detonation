#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from corner_radius_study import radius_tag
from corner_solution_frames import pressure_threshold_from_csv


def parse_radii(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def load_final(case_dir: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    data = np.loadtxt(case_dir / "final_results.csv", delimiter=",", skiprows=1, usecols=(0, 1, 5, 6, 12))
    return data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4] != 0.0


def grid_from_points(x: np.ndarray, y: np.ndarray, values: np.ndarray, solid: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_levels = np.unique(x)
    y_levels = np.unique(y)
    dx = x_levels[1] - x_levels[0] if len(x_levels) > 1 else 1.0
    dy = y_levels[1] - y_levels[0] if len(y_levels) > 1 else 1.0
    ix = np.rint((x - x_levels[0]) / dx).astype(int)
    iy = np.rint((y - y_levels[0]) / dy).astype(int)
    grid = np.full((len(y_levels), len(x_levels)), np.nan)
    fluid = ~solid
    grid[iy[fluid], ix[fluid]] = values[fluid]
    return x_levels, y_levels, grid


def dead_mask_from_summary(case_dir: Path, x: np.ndarray, y: np.ndarray, w: np.ndarray, solid: np.ndarray) -> np.ndarray:
    p_threshold = pressure_threshold_from_csv(case_dir / "step_0_initial.csv", 0.01)
    max_pressure = np.full_like(w, -math.inf)
    frames = [case_dir / "step_0_initial.csv"] + sorted(case_dir.glob("step_*_time_*.csv")) + [case_dir / "final_results.csv"]
    for frame in frames:
        data = np.loadtxt(frame, delimiter=",", skiprows=1, usecols=(5, 12))
        fluid = data[:, 1] == 0.0
        max_pressure[fluid] = np.maximum(max_pressure[fluid], data[fluid, 0])
    analysis = (~solid) & (x >= 2.0) & (y <= 3.0) & (max_pressure >= p_threshold)
    return analysis & (w > 0.05)


def plot_field_collection(output_root: Path, radii: list[float], out_dir: Path, field: str) -> None:
    loaded = []
    all_pressure = []
    for radius in radii:
        case_dir = output_root / radius_tag(radius)
        x, y, p, w, solid = load_final(case_dir)
        values = p if field == "p" else w
        if field == "p":
            all_pressure.extend(values[~solid & np.isfinite(values)].tolist())
        loaded.append((radius, x, y, values, solid))

    if field == "p":
        finite = sorted(value for value in all_pressure if math.isfinite(value))
        vmin = 0.0
        vmax = finite[int(0.995 * (len(finite) - 1))] if finite else None
        cmap_name = "inferno"
        label = "P [Mbar]"
        title = "Final pressure"
        filename = "final_pressure_by_radius.png"
    else:
        vmin = 0.0
        vmax = 1.0
        cmap_name = "viridis_r"
        label = "W"
        title = "Final W"
        filename = "final_w_by_radius.png"

    fig, axes = plt.subplots(1, len(radii), figsize=(3.4 * len(radii), 3.0), sharex=True, sharey=True, constrained_layout=True)
    if len(radii) == 1:
        axes = [axes]
    image = None
    for axis, (radius, x, y, values, solid) in zip(axes, loaded):
        x_levels, y_levels, grid = grid_from_points(x, y, values, solid)
        extent = [float(x_levels.min()), float(x_levels.max()), float(y_levels.min()), float(y_levels.max())]
        cmap = plt.get_cmap(cmap_name).copy()
        cmap.set_bad("#d9d9d9")
        image = axis.imshow(grid, extent=extent, origin="lower", aspect="equal", cmap=cmap, vmin=vmin, vmax=vmax)
        axis.set_title(f"R={radius:g} cm", fontsize=10)
        axis.set_xlim(0.0, 7.0)
        axis.set_ylim(0.0, 4.0)
        axis.set_xlabel("x [cm]")
    axes[0].set_ylabel("y [cm]")
    if image is not None:
        fig.colorbar(image, ax=axes, shrink=0.82, pad=0.012, label=label)
    fig.suptitle(title, fontsize=12)
    fig.savefig(out_dir / filename, dpi=190)
    plt.close(fig)


def plot_dead_collection(output_root: Path, radii: list[float], out_dir: Path) -> None:
    fig, axes = plt.subplots(1, len(radii), figsize=(3.4 * len(radii), 3.0), sharex=True, sharey=True, constrained_layout=True)
    if len(radii) == 1:
        axes = [axes]
    image = None
    for axis, radius in zip(axes, radii):
        case_dir = output_root / radius_tag(radius)
        x, y, _, w, solid = load_final(case_dir)
        dead = dead_mask_from_summary(case_dir, x, y, w, solid)
        values = np.full_like(w, np.nan)
        values[(~solid) & (x >= 2.0) & (y <= 3.0)] = 0.0
        values[dead] = 1.0
        x_levels, y_levels, grid = grid_from_points(x, y, values, solid)
        extent = [float(x_levels.min()), float(x_levels.max()), float(y_levels.min()), float(y_levels.max())]
        cmap = plt.get_cmap("gray_r").copy()
        cmap.set_bad("#d9d9d9")
        image = axis.imshow(grid, extent=extent, origin="lower", aspect="equal", cmap=cmap, vmin=0.0, vmax=1.0)
        axis.set_title(f"R={radius:g} cm", fontsize=10)
        axis.set_xlim(0.0, 7.0)
        axis.set_ylim(0.0, 4.0)
        axis.set_xlabel("x [cm]")
    axes[0].set_ylabel("y [cm]")
    if image is not None:
        fig.colorbar(image, ax=axes, shrink=0.82, pad=0.012, label="dead = 1")
    fig.suptitle("Binary dead-zone map", fontsize=12)
    fig.savefig(out_dir / "dead_zone_by_radius.png", dpi=190)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Regenerate report field figures with equal physical aspect ratio.")
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--radii", default="0,0.75,1,2")
    parser.add_argument("--out-dir", default="tests/mader_2de/report/figures")
    args = parser.parse_args()

    output_root = Path(args.output_root).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    radii = parse_radii(args.radii)
    plot_field_collection(output_root, radii, out_dir, "p")
    plot_field_collection(output_root, radii, out_dir, "w")
    plot_dead_collection(output_root, radii, out_dir)
    print(f"[write] {out_dir / 'final_pressure_by_radius.png'}")
    print(f"[write] {out_dir / 'final_w_by_radius.png'}")
    print(f"[write] {out_dir / 'dead_zone_by_radius.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
