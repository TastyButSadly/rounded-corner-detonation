#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import io
import math
from pathlib import Path

import imageio.v2 as imageio
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


def f(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return math.nan


def load_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def physical_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in rows if int(float(row["is_solid"])) == 0]


def frame_time(path: Path) -> float:
    stem = path.stem
    if stem == "step_0_initial":
        return 0.0
    if stem == "final_results":
        return math.inf
    for marker in ("_time_", "_t_"):
        if marker in stem:
            try:
                return float(stem.rsplit(marker, 1)[1])
            except Exception:
                return math.nan
    return math.nan


def sorted_frames(case_dir: Path) -> list[Path]:
    step_frames = sorted(
        list(case_dir.glob("step_*_time_*.csv")) + list(case_dir.glob("step_*_t_*.csv")),
        key=frame_time,
    )
    frames: list[Path] = []
    initial = case_dir / "step_0_initial.csv"
    final = case_dir / "final_results.csv"
    if initial.exists():
        frames.append(initial)
    frames.extend(step_frames)
    if final.exists():
        frames.append(final)
    return frames


def build_grid(rows: list[dict[str, str]], field: str) -> tuple[list[float], list[float], np.ndarray]:
    x_levels = sorted({f(row["x"]) for row in rows})
    y_levels = sorted({f(row["y"]) for row in rows})
    x_index = {value: idx for idx, value in enumerate(x_levels)}
    y_index = {value: idx for idx, value in enumerate(y_levels)}
    grid = np.full((len(y_levels), len(x_levels)), np.nan, dtype=float)
    for row in rows:
        x = f(row["x"])
        y = f(row["y"])
        value = f(row[field])
        if math.isfinite(x) and math.isfinite(y) and math.isfinite(value):
            grid[y_index[y], x_index[x]] = value
    return x_levels, y_levels, grid


def render_frame(
    rows: list[dict[str, str]],
    *,
    time_value: float,
    pressure_vmax: float,
    pressure_cmap: str,
    mode: str,
) -> Image.Image:
    x_levels, y_levels, p_grid = build_grid(rows, "p")
    _, _, w_grid = build_grid(rows, "w")
    extent = [min(x_levels), max(x_levels), min(y_levels), max(y_levels)]

    if mode == "pressure":
        fig, axes = plt.subplots(1, 1, figsize=(7.5, 4.8), sharex=True, sharey=True)
        axes = [axes]
        images = [
            (axes[0], p_grid, "Pressure P [Mbar]", pressure_cmap, 0.0, pressure_vmax),
        ]
    elif mode == "w":
        fig, axes = plt.subplots(1, 1, figsize=(7.5, 4.8), sharex=True, sharey=True)
        axes = [axes]
        images = [
            (axes[0], w_grid, "Unreacted mass fraction W [-]", "viridis_r", 0.0, 1.0),
        ]
    else:
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), sharex=True, sharey=True)
        images = [
            (axes[0], p_grid, "Pressure P [Mbar]", pressure_cmap, 0.0, pressure_vmax),
            (axes[1], w_grid, "Unreacted mass fraction W [-]", "viridis_r", 0.0, 1.0),
        ]

    for axis, grid, title, cmap, vmin, vmax in images:
        image = axis.imshow(
            grid,
            extent=extent,
            origin="lower",
            aspect="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        axis.set_title(title)
        axis.set_xlabel("x [cm]")
        fig.colorbar(image, ax=axis, shrink=0.9)
    axes[0].set_ylabel("y [cm]")

    if math.isfinite(time_value):
        title = f"Corner turning evolution, t = {time_value:.3f} μs"
    else:
        title = "Corner turning evolution, final exported state"
    fig.suptitle(title, fontsize=14)
    fig.text(
        0.5,
        0.01,
        "Units: x, y [cm]; t [μs]; ρ [g/cm³]; P [Mbar]; U, V [cm/μs]; T [K]; W [-]",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    fig.tight_layout(rect=(0.03, 0.05, 0.98, 0.94))

    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=130)
    plt.close(fig)
    buffer.seek(0)
    return Image.open(buffer).convert("RGB")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build a GIF for corner-turning outputs")
    parser.add_argument("--output-root", required=True, help="Directory that contains the corner case directory")
    parser.add_argument("--case", default="corner", help="Case subdirectory name inside output-root")
    parser.add_argument("--stride", type=int, default=1, help="Use every N-th saved frame")
    parser.add_argument("--duration-ms", type=int, default=140, help="Per-frame duration in milliseconds")
    parser.add_argument("--mode", choices=["both", "pressure", "w"], default="both")
    parser.add_argument("--pressure-percentile", type=float, default=100.0, help="Percentile-based cap for pressure color scale")
    parser.add_argument("--pressure-vmax", type=float, default=None, help="Explicit upper bound for pressure color scale")
    parser.add_argument("--pressure-cmap", type=str, default="turbo", help="Matplotlib colormap for pressure panels")
    args = parser.parse_args()

    case_dir = Path(args.output_root).resolve() / args.case
    if not case_dir.exists():
        raise SystemExit(f"Missing case directory: {case_dir}")

    frames = sorted_frames(case_dir)
    if not frames:
        raise SystemExit(f"No frames found in {case_dir}")
    stride = max(args.stride, 1)
    selected = [frame for index, frame in enumerate(frames) if index % stride == 0]
    if selected[-1] != frames[-1]:
        selected.append(frames[-1])

    pressure_samples: list[float] = []
    for csv_path in selected:
        for row in physical_rows(load_rows(csv_path)):
            value = f(row["p"])
            if math.isfinite(value):
                pressure_samples.append(value)
    if args.pressure_vmax is not None:
        pressure_vmax = args.pressure_vmax
    elif pressure_samples:
        percentile = min(max(args.pressure_percentile, 0.0), 100.0)
        pressure_vmax = float(np.percentile(np.asarray(pressure_samples), percentile))
    else:
        pressure_vmax = 1.0
    pressure_vmax = max(pressure_vmax, 1e-6)

    images = []
    for csv_path in selected:
        rows = physical_rows(load_rows(csv_path))
        images.append(
            np.asarray(
                render_frame(
                    rows,
                    time_value=frame_time(csv_path),
                    pressure_vmax=pressure_vmax,
                    pressure_cmap=args.pressure_cmap,
                    mode=args.mode,
                )
            )
        )

    plot_dir = case_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    suffix = {
        "both": "corner_evolution.gif",
        "pressure": "corner_pressure.gif",
        "w": "corner_w.gif",
    }[args.mode]
    gif_path = plot_dir / suffix
    imageio.mimsave(gif_path, images, duration=max(args.duration_ms, 20) / 1000.0, loop=0)
    print(gif_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
