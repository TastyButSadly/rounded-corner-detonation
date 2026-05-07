#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

from corner_radius_study import f, frame_time, is_solid, load_rows, radius_tag, sorted_frames


FIELD_COLUMNS = {
    "p": 5,
    "w": 6,
}
OVERLAY_DILATION_CELLS = 3


def parse_values(raw: str, cast=float) -> list:
    values = []
    for item in raw.split(","):
        item = item.strip()
        if item:
            values.append(cast(item))
    return values


def nearest_frame(frames: list[Path], target: float) -> Path:
    finite = [path for path in frames if math.isfinite(frame_time(path))]
    if not finite:
        raise SystemExit("No finite-time CSV frames found.")
    return min(finite, key=lambda path: abs(frame_time(path) - target))


def infer_levels(rows: list[dict[str, str]]) -> tuple[list[float], list[float], dict[tuple[str, str], tuple[int, int]]]:
    x_levels = sorted({f(row["x"]) for row in rows if math.isfinite(f(row["x"]))})
    y_levels = sorted({f(row["y"]) for row in rows if math.isfinite(f(row["y"]))})
    x_index = {value: idx for idx, value in enumerate(x_levels)}
    y_index = {value: idx for idx, value in enumerate(y_levels)}
    index = {}
    for row in rows:
        x = f(row["x"])
        y = f(row["y"])
        if math.isfinite(x) and math.isfinite(y):
            index[(row["x"], row["y"])] = (y_index[y], x_index[x])
    return x_levels, y_levels, index


def field_grid(
    rows: list[dict[str, str]],
    field: str,
    x_levels: list[float],
    y_levels: list[float],
    index: dict[tuple[str, str], tuple[int, int]],
) -> list[list[float]]:
    grid = [[math.nan for _ in x_levels] for _ in y_levels]
    for row in rows:
        key = (row["x"], row["y"])
        if key not in index or is_solid(row):
            continue
        value = f(row[field])
        if math.isfinite(value):
            iy, ix = index[key]
            grid[iy][ix] = value
    return grid


def dead_overlay_grid(
    rows: list[dict[str, str]],
    arrived_keys: set[tuple[str, str]],
    w_threshold: float,
    x_analysis_min: float,
    y_analysis_max: float,
    x_levels: list[float],
    y_levels: list[float],
    index: dict[tuple[str, str], tuple[int, int]],
) -> list[list[float]]:
    grid = [[math.nan for _ in x_levels] for _ in y_levels]
    for row in rows:
        key = (row["x"], row["y"])
        if key not in index or key not in arrived_keys or is_solid(row):
            continue
        x = f(row["x"])
        y = f(row["y"])
        if x + 1e-12 < x_analysis_min or y - 1e-12 > y_analysis_max:
            continue
        if f(row["w"]) > w_threshold:
            iy, ix = index[key]
            grid[iy][ix] = 1.0
    return grid


def pressure_reference(initial_rows: list[dict[str, str]], arrival_fraction: float) -> tuple[float, float, float]:
    values = [f(row["p"]) for row in initial_rows if not is_solid(row) and math.isfinite(f(row["p"]))]
    if not values:
        raise SystemExit("Could not infer pressure reference from initial frame.")
    p0 = min(values)
    p_ref = max(values)
    return p0, p_ref, p0 + arrival_fraction * (p_ref - p0)


def pressure_threshold_from_csv(csv_path: Path, arrival_fraction: float) -> float:
    data = np.loadtxt(csv_path, delimiter=",", skiprows=1, usecols=(5, 12))
    fluid = data[:, 1] == 0.0
    values = data[fluid, 0]
    p0 = float(np.min(values))
    p_ref = float(np.max(values))
    return p0 + arrival_fraction * (p_ref - p0)


def case_frame_data(
    frames: list[Path],
    target_frames: list[Path],
    p_threshold: float,
    use_arrival_history: bool,
) -> list[tuple[Path, list[dict[str, str]], set[tuple[str, str]]]]:
    if not use_arrival_history:
        result = []
        for path in target_frames:
            rows = load_rows(path)
            visible_keys = {(row["x"], row["y"]) for row in rows if not is_solid(row)}
            result.append((path, rows, visible_keys))
        return result

    targets = {path.resolve() for path in target_frames}
    max_time = max(frame_time(path) for path in target_frames)

    arrived: set[tuple[str, str]] = set()
    result: dict[Path, tuple[Path, list[dict[str, str]], set[tuple[str, str]]]] = {}
    for path in frames:
        time = frame_time(path)
        if not math.isfinite(time) or time > max_time + 1e-12:
            continue
        rows = load_rows(path)
        for row in rows:
            if is_solid(row):
                continue
            if f(row["p"]) >= p_threshold:
                arrived.add((row["x"], row["y"]))
        resolved = path.resolve()
        if resolved in targets:
            result[resolved] = (path, rows, set(arrived))
    return [result[path.resolve()] for path in target_frames]


def load_frame_arrays(
    csv_path: Path,
    field: str,
    w_threshold: float,
    p_threshold: float,
    x_analysis_min: float,
    y_analysis_max: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    value_col = FIELD_COLUMNS[field]
    if field == "w":
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1, usecols=(0, 1, value_col, 5, 12))
        x = data[:, 0]
        y = data[:, 1]
        values = data[:, 2]
        w_values = values
        p_values = data[:, 3]
        solid = data[:, 4] != 0.0
    else:
        data = np.loadtxt(csv_path, delimiter=",", skiprows=1, usecols=(0, 1, value_col, 6, 12))
        x = data[:, 0]
        y = data[:, 1]
        values = data[:, 2]
        w_values = data[:, 3]
        p_values = values
        solid = data[:, 4] != 0.0

    x_levels = np.unique(x)
    y_levels = np.unique(y)
    dx = x_levels[1] - x_levels[0] if len(x_levels) > 1 else 1.0
    dy = y_levels[1] - y_levels[0] if len(y_levels) > 1 else 1.0
    ix = np.rint((x - x_levels[0]) / dx).astype(int)
    iy = np.rint((y - y_levels[0]) / dy).astype(int)

    grid = np.full((len(y_levels), len(x_levels)), np.nan)
    overlay = np.full_like(grid, np.nan)
    fluid = ~solid
    grid[iy[fluid], ix[fluid]] = values[fluid]
    dead = fluid & (x >= x_analysis_min) & (y <= y_analysis_max) & (p_values >= p_threshold) & (w_values > w_threshold)
    dead_grid = np.zeros_like(grid, dtype=bool)
    dead_grid[iy[dead], ix[dead]] = True
    if OVERLAY_DILATION_CELLS > 0:
        base = dead_grid.copy()
        for y_shift in range(-OVERLAY_DILATION_CELLS, OVERLAY_DILATION_CELLS + 1):
            for x_shift in range(-OVERLAY_DILATION_CELLS, OVERLAY_DILATION_CELLS + 1):
                if x_shift == 0 and y_shift == 0:
                    continue
                shifted = np.roll(np.roll(base, y_shift, axis=0), x_shift, axis=1)
                if y_shift > 0:
                    shifted[:y_shift, :] = False
                elif y_shift < 0:
                    shifted[y_shift:, :] = False
                if x_shift > 0:
                    shifted[:, :x_shift] = False
                elif x_shift < 0:
                    shifted[:, x_shift:] = False
                dead_grid |= shifted
    overlay[dead_grid] = 1.0
    return x_levels, y_levels, grid, overlay


def plot_montage(
    output_root: Path,
    radii: list[float],
    times: list[float],
    field: str,
    out_path: Path,
    w_threshold: float,
    arrival_fraction: float,
    x_analysis_min: float,
    y_analysis_max: float,
    use_arrival_history: bool,
    show_overlay: bool,
) -> None:
    cases = []
    all_values = []
    for radius in radii:
        case_dir = output_root / radius_tag(radius)
        frames = sorted_frames(case_dir)
        if not frames:
            raise SystemExit(f"No frames found for {case_dir}")
        p_threshold = pressure_threshold_from_csv(case_dir / "step_0_initial.csv", arrival_fraction)
        target_frames = [nearest_frame(frames, target) for target in times]
        case_frames = []
        for path in target_frames:
            x_levels, y_levels, grid, overlay = load_frame_arrays(
                path,
                field,
                w_threshold,
                p_threshold,
                x_analysis_min,
                y_analysis_max,
            )
            case_frames.append((path, x_levels, y_levels, grid, overlay))
        if field == "p":
            for _, _, _, grid, _ in case_frames:
                all_values.extend(grid[np.isfinite(grid)].tolist())
        cases.append((radius, case_frames))

    if field == "p":
        finite = sorted(value for value in all_values if math.isfinite(value))
        vmin = 0.0
        vmax = finite[int(0.995 * (len(finite) - 1))] if finite else None
        cmap_name = "inferno"
        colorbar_label = "P [Mbar]"
    else:
        vmin = 0.0
        vmax = 1.0
        cmap_name = "viridis_r"
        colorbar_label = "W"

    fig, axes = plt.subplots(
        len(radii),
        len(times),
        figsize=(3.05 * len(times), 2.05 * len(radii)),
        sharex=True,
        sharey=True,
    )
    if len(radii) == 1:
        axes = [axes]
    if len(times) == 1:
        axes = [[axis] for axis in axes]

    dead_cmap = ListedColormap(["#d62728"])
    image = None
    for row_axes, (radius, case_frames) in zip(axes, cases):
        for axis, (path, x_levels, y_levels, grid, overlay) in zip(row_axes, case_frames):
            extent = [min(x_levels), max(x_levels), min(y_levels), max(y_levels)]
            cmap = plt.get_cmap(cmap_name).copy()
            cmap.set_bad("#d8d8d8")
            image = axis.imshow(
                grid,
                extent=extent,
                origin="lower",
                aspect="equal",
                cmap=cmap,
                vmin=vmin,
                vmax=vmax,
            )
            if show_overlay:
                axis.imshow(overlay, extent=extent, origin="lower", aspect="equal", cmap=dead_cmap, alpha=0.58)
            axis.set_title(f"R={radius:g}, t={frame_time(path):.2f} us", fontsize=9)
            axis.set_xlim(1.6, 6.8)
            axis.set_ylim(0.0, 4.0)
            axis.tick_params(labelsize=8)

    for axis in axes[-1]:
        axis.set_xlabel("x [cm]", fontsize=9)
    for row_axes in axes:
        row_axes[0].set_ylabel("y [cm]", fontsize=9)
    if image is not None:
        cbar = fig.colorbar(image, ax=[axis for row_axes in axes for axis in row_axes], shrink=0.86, pad=0.012)
        cbar.set_label(colorbar_label)
    title = f"{field.upper()} evolution"
    if show_overlay:
        title += " with dead-zone candidates in red"
    fig.suptitle(title, y=0.995, fontsize=12)
    fig.savefig(out_path, dpi=190, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Build solution-frame montages with dead-zone overlays.")
    parser.add_argument("--output-root", required=True)
    parser.add_argument("--radii", default="0,0.75,1,2")
    parser.add_argument("--times", default="1.5,2.25,3.0,3.75,4.5,5.5")
    parser.add_argument("--out-dir", default="tests/mader_2de/report/figures")
    parser.add_argument("--w-threshold", type=float, default=0.05)
    parser.add_argument("--arrival-fraction", type=float, default=0.01)
    parser.add_argument("--arrival-history-overlay", action="store_true", help="Use cumulative pressure-arrival history for the red overlay. Slower for fine grids.")
    parser.add_argument("--show-overlay", action="store_true", help="Draw red W-threshold overlay on every time frame.")
    parser.add_argument("--x-analysis-min", type=float, default=2.0)
    parser.add_argument("--y-analysis-max", type=float, default=3.0)
    args = parser.parse_args()

    output_root = Path(args.output_root).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    radii = parse_values(args.radii, float)
    times = parse_values(args.times, float)

    label = output_root.name.replace("output_corner_radius_convergence_", "")
    plot_montage(
        output_root,
        radii,
        times,
        "p",
        out_dir / f"solution_evolution_pressure_dead_{label}.png",
        args.w_threshold,
        args.arrival_fraction,
        args.x_analysis_min,
        args.y_analysis_max,
        args.arrival_history_overlay,
        args.show_overlay,
    )
    plot_montage(
        output_root,
        radii,
        times,
        "w",
        out_dir / f"solution_evolution_w_dead_{label}.png",
        args.w_threshold,
        args.arrival_fraction,
        args.x_analysis_min,
        args.y_analysis_max,
        args.arrival_history_overlay,
        args.show_overlay,
    )
    print(f"[write] {out_dir / f'solution_evolution_pressure_dead_{label}.png'}")
    print(f"[write] {out_dir / f'solution_evolution_w_dead_{label}.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
