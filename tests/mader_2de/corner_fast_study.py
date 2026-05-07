#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from corner_radius_study import (
    DEFAULT_RADII,
    format_number,
    generate_config,
    parse_radii,
    radius_tag,
    resolve_path,
    resolve_solver,
    result_is_critical,
    run_solver,
    sorted_frames,
)


@dataclass(frozen=True)
class FastResult:
    radius: float
    s_dead: float
    s_hat: float
    s_analysis: float
    s_box: float
    arrived_fraction: float
    arrived_cells: int
    dead_cells: int


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


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


def write_field_plot(path: Path, x: np.ndarray, y: np.ndarray, values: np.ndarray, solid: np.ndarray, title: str, cmap_name: str, vmin: float | None, vmax: float | None) -> None:
    x_levels, y_levels, grid = grid_from_points(x, y, values, solid)
    extent = [float(x_levels.min()), float(x_levels.max()), float(y_levels.min()), float(y_levels.max())]
    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad("#d9d9d9")
    fig, axis = plt.subplots(figsize=(7.8, 4.8))
    image = axis.imshow(grid, extent=extent, origin="lower", aspect="equal", cmap=cmap, vmin=vmin, vmax=vmax)
    axis.set_title(title)
    axis.set_xlabel("x [cm]")
    axis.set_ylabel("y [cm]")
    fig.colorbar(image, ax=axis, shrink=0.82)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def write_dead_plot(path: Path, x: np.ndarray, y: np.ndarray, solid: np.ndarray, analysis: np.ndarray, dead: np.ndarray, title: str) -> None:
    values = np.full_like(x, np.nan, dtype=float)
    values[analysis] = 0.0
    values[dead] = 1.0
    x_levels, y_levels, grid = grid_from_points(x, y, values, solid & ~analysis)
    extent = [float(x_levels.min()), float(x_levels.max()), float(y_levels.min()), float(y_levels.max())]
    cmap = plt.get_cmap("gray_r").copy()
    cmap.set_bad("#d9d9d9")
    fig, axis = plt.subplots(figsize=(7.8, 4.8))
    image = axis.imshow(grid, extent=extent, origin="lower", aspect="equal", cmap=cmap, vmin=0.0, vmax=1.0)
    axis.set_title(title)
    axis.set_xlabel("x [cm]")
    axis.set_ylabel("y [cm]")
    fig.colorbar(image, ax=axis, shrink=0.82, label="dead = 1")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def analyze_case_fast(
    case_dir: Path,
    radius: float,
    x_analysis_min: float,
    y_analysis_max: float,
    w_threshold: float,
) -> FastResult:
    initial = case_dir / "step_0_initial.csv"
    final = case_dir / "final_results.csv"
    if not initial.exists() or not final.exists():
        raise SystemExit(f"Missing initial or final CSV in {case_dir}")

    initial_data = np.loadtxt(initial, delimiter=",", skiprows=1, usecols=(5, 12))
    initial_fluid = initial_data[:, 1] == 0.0
    p0 = float(np.min(initial_data[initial_fluid, 0]))
    p_ref = float(np.max(initial_data[initial_fluid, 0]))
    p_threshold = p0 + 0.01 * (p_ref - p0)

    final_data = np.loadtxt(final, delimiter=",", skiprows=1, usecols=(0, 1, 5, 6, 12))
    x = final_data[:, 0]
    y = final_data[:, 1]
    p_final = final_data[:, 2]
    w_final = final_data[:, 3]
    solid = final_data[:, 4] != 0.0
    x_levels = np.unique(x)
    y_levels = np.unique(y)
    dx = float(np.min(np.diff(x_levels)))
    dy = float(np.min(np.diff(y_levels)))
    cell_area = dx * dy

    max_pressure = np.full_like(p_final, -math.inf, dtype=float)
    for frame in sorted_frames(case_dir):
        data = np.loadtxt(frame, delimiter=",", skiprows=1, usecols=(5, 12))
        fluid = data[:, 1] == 0.0
        max_pressure[fluid] = np.maximum(max_pressure[fluid], data[fluid, 0])

    box = (~solid) & (x >= x_analysis_min - 1e-12) & (y <= y_analysis_max + 1e-12)
    analysis = box & (max_pressure >= p_threshold)
    dead = analysis & (w_final > w_threshold)

    box_cells = int(np.count_nonzero(box))
    arrived_cells = int(np.count_nonzero(analysis))
    dead_cells = int(np.count_nonzero(dead))
    s_box = box_cells * cell_area
    s_analysis = arrived_cells * cell_area
    s_dead = dead_cells * cell_area
    s_hat = s_dead / s_analysis if s_analysis > 0.0 else math.nan
    arrived_fraction = arrived_cells / box_cells if box_cells > 0 else math.nan

    plot_dir = case_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    tag = radius_tag(radius)
    write_field_plot(plot_dir / f"{tag}_final_pressure.png", x, y, p_final, solid, f"Final pressure, R = {radius:g} cm", "inferno", None, None)
    write_field_plot(plot_dir / f"{tag}_final_w.png", x, y, w_final, solid, f"Final W, R = {radius:g} cm", "viridis_r", 0.0, 1.0)
    write_dead_plot(plot_dir / f"{tag}_dead_zone.png", x, y, solid, analysis, dead, f"Dead-zone map, R = {radius:g} cm")
    return FastResult(radius, s_dead, s_hat, s_analysis, s_box, arrived_fraction, arrived_cells, dead_cells)


def write_summary(output_root: Path, results: list[FastResult], epsilon: float) -> None:
    summary_path = output_root / "summary.csv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "R_corner_cm",
            "S_dead_cm2",
            "S_hat",
            "S_analysis_cm2",
            "S_box_cm2",
            "arrived_fraction",
            "arrived_cells",
            "dead_cells",
        ])
        for result in results:
            writer.writerow([
                format_number(result.radius),
                format_number(result.s_dead),
                format_number(result.s_hat),
                format_number(result.s_analysis),
                format_number(result.s_box),
                format_number(result.arrived_fraction),
                result.arrived_cells,
                result.dead_cells,
            ])

    critical = next((result.radius for result in sorted(results, key=lambda item: item.radius) if result_is_critical(result, epsilon)), None)
    if critical is None:
        text = f"R_cr > {max(result.radius for result in results):g} cm for epsilon = {epsilon:g}\n"
    else:
        text = f"R_cr = {critical:g} cm for epsilon = {epsilon:g}\n"
    (output_root / "critical_radius.txt").write_text(text, encoding="utf-8")


def plot_summary(output_root: Path, results: list[FastResult]) -> None:
    sorted_results = sorted(results, key=lambda item: item.radius)
    radii = [result.radius for result in sorted_results]
    s_hat = [result.s_hat for result in sorted_results]
    s_dead = [result.s_dead for result in sorted_results]
    plot_dir = output_root / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    for values, ylabel, filename in (
        (s_dead, "S_dead [cm^2]", "S_dead_vs_R.png"),
        (s_hat, "S_dead / S_analysis [-]", "S_hat_vs_R.png"),
    ):
        fig, axis = plt.subplots(figsize=(7.5, 4.8))
        axis.plot(radii, values, marker="o", linewidth=1.8)
        axis.set_xlabel("R_corner [cm]")
        axis.set_ylabel(ylabel)
        axis.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(plot_dir / filename, dpi=170)
        plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run and analyze fine rounded-corner studies using numpy-based post-processing.")
    parser.add_argument("--solver", default=None)
    parser.add_argument("--base-config", default="tests/mader_2de/configs/corner_radius_study_base.ini")
    parser.add_argument("--output-root", default="tests/mader_2de/output_corner_radius_fine")
    parser.add_argument("--radii", default=DEFAULT_RADII)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--tmax", type=float, default=5.5)
    parser.add_argument("--dt-out", type=float, default=0.25)
    parser.add_argument("--x-analysis-min", type=float, default=2.0)
    parser.add_argument("--y-analysis-max", type=float, default=3.0)
    parser.add_argument("--w-threshold", type=float, default=0.05)
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--skip-run", action="store_true")
    parser.add_argument("--keep-existing", action="store_true")
    args = parser.parse_args()

    root = repo_root()
    output_root = resolve_path(root, args.output_root)
    base_config = resolve_path(root, args.base_config)
    output_root.mkdir(parents=True, exist_ok=True)
    solver = None if args.skip_run else resolve_solver(root, args.solver)
    radii = parse_radii(args.radii)
    results: list[FastResult] = []
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str((root / "tests" / "mader_2de" / ".mplcache").resolve()))

    for radius in radii:
        tag = radius_tag(radius)
        config_path = output_root / "configs" / f"corner_{tag}.ini"
        case_dir = output_root / tag
        generate_config(base_config, config_path, radius, tag, args.nx, args.ny, args.tmax, args.dt_out)
        if solver is not None:
            run_solver(solver, config_path, case_dir, root, args.keep_existing)
        result = analyze_case_fast(case_dir, radius, args.x_analysis_min, args.y_analysis_max, args.w_threshold)
        results.append(result)
        print(
            f"[analyze] R={radius:g}: S_dead={result.s_dead:.6g}, "
            f"S_hat={result.s_hat:.6g}, dead_cells={result.dead_cells}, arrived_cells={result.arrived_cells}",
            flush=True,
        )

    write_summary(output_root, results, args.epsilon)
    plot_summary(output_root, results)
    print(f"[summary] wrote {output_root / 'summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
