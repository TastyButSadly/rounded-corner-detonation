#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


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


def by_y(rows: list[dict[str, str]]) -> dict[float, list[dict[str, str]]]:
    grouped: dict[float, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(f(row["y"]), []).append(row)
    return grouped


def pick_row_by_y(rows: list[dict[str, str]], target_y: float | None = None) -> list[dict[str, str]]:
    grouped = by_y(rows)
    if not grouped:
        return []
    y_levels = sorted(grouped)
    if target_y is None:
        target_y = y_levels[len(y_levels) // 2]
    best_y = min(y_levels, key=lambda y: abs(y - target_y))
    return sorted(grouped[best_y], key=lambda row: f(row["x"]))


def frame_time(csv_path: Path) -> float:
    name = csv_path.name
    if name == "step_0_initial.csv":
        return 0.0
    if name == "final_results.csv":
        return math.inf
    for marker in ("_time_", "_t_"):
        if marker in name:
            try:
                return float(name.split(marker, 1)[1].rsplit(".csv", 1)[0])
            except ValueError:
                return math.nan
    return math.nan


def sorted_frames(case_dir: Path) -> list[Path]:
    step_frames = sorted(
        list(case_dir.glob("step_*_time_*.csv")) + list(case_dir.glob("step_*_t_*.csv")),
        key=frame_time,
    )
    initial = case_dir / "step_0_initial.csv"
    final = case_dir / "final_results.csv"
    frames: list[Path] = []
    if initial.exists():
        frames.append(initial)
    frames.extend(step_frames)
    if final.exists():
        frames.append(final)
    return frames


def front_position_from_w(profile: list[dict[str, str]], threshold: float = 0.5) -> float:
    front = math.nan
    for row in profile:
        x = f(row["x"])
        w = f(row["w"])
        if math.isfinite(x) and math.isfinite(w) and w <= threshold:
            front = x
        else:
            break
    return front


def ensure_plot_dir(case_dir: Path) -> Path:
    plot_dir = case_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    return plot_dir


def annotate_units(fig: plt.Figure, extra: str = "") -> None:
    units = "Units: x, y [cm]; t [μs]; ρ [g/cm³]; P [Mbar]; U, V [cm/μs]; T [K]; W [-]"
    text = units if not extra else f"{units}\n{extra}"
    fig.text(0.5, 0.01, text, ha="center", va="bottom", fontsize=9)


def plot_sod(case_dir: Path) -> None:
    final_rows = physical_rows(load_rows(case_dir / "final_results.csv"))
    profile = pick_row_by_y(final_rows)
    if not profile:
        return

    x = [f(row["x"]) for row in profile]
    rho = [f(row["rho"]) for row in profile]
    u = [f(row["u"]) for row in profile]
    p = [f(row["p"]) for row in profile]
    rho_exact = [f(row["rho_exact"]) for row in profile]
    u_exact = [f(row["u_exact"]) for row in profile]
    p_exact = [f(row["p_exact"]) for row in profile]

    fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
    series = [
        (axes[0], rho, rho_exact, "Density ρ [g/cm³]"),
        (axes[1], p, p_exact, "Pressure P [Mbar]"),
        (axes[2], u, u_exact, "Velocity U [cm/μs]"),
    ]
    for axis, numeric, exact, ylabel in series:
        axis.plot(x, numeric, linewidth=2, label="2DE result")
        if any(math.isfinite(value) for value in exact):
            axis.plot(x, exact, "--", linewidth=1.5, label="Reference")
        axis.set_ylabel(ylabel)
        axis.grid(True, alpha=0.3)
        axis.legend(loc="best")

    axes[-1].set_xlabel("x [cm]")
    fig.suptitle("Sod Shock Tube, central row", fontsize=14)
    annotate_units(fig, "Reference curves are taken from the exact solution columns exported by the solver.")
    fig.tight_layout(rect=(0.04, 0.05, 0.98, 0.97))
    fig.savefig(ensure_plot_dir(case_dir) / "sod_profiles.png", dpi=160)
    plt.close(fig)


def plot_cj(case_dir: Path) -> None:
    final_rows = physical_rows(load_rows(case_dir / "final_results.csv"))
    profile = pick_row_by_y(final_rows)
    if not profile:
        return

    x = [f(row["x"]) for row in profile]
    rho = [f(row["rho"]) for row in profile]
    u = [f(row["u"]) for row in profile]
    p = [f(row["p"]) for row in profile]
    w = [f(row["w"]) for row in profile]
    temperature = [f(row["temperature"]) for row in profile]

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    panels = [
        (axes[0][0], p, "Pressure P [Mbar]"),
        (axes[0][1], rho, "Density ρ [g/cm³]"),
        (axes[1][0], w, "Unreacted mass fraction W [-]"),
        (axes[1][1], temperature, "Temperature T [K]"),
    ]
    for axis, values, ylabel in panels:
        axis.plot(x, values, linewidth=2)
        axis.set_ylabel(ylabel)
        axis.grid(True, alpha=0.3)
    axes[1][0].set_xlabel("x [cm]")
    axes[1][1].set_xlabel("x [cm]")
    fig.suptitle("CJ-style detonation, central row", fontsize=14)
    annotate_units(fig)
    fig.tight_layout(rect=(0.04, 0.05, 0.98, 0.95))
    plot_dir = ensure_plot_dir(case_dir)
    fig.savefig(plot_dir / "cj_profiles.png", dpi=160)
    plt.close(fig)

    frames = sorted_frames(case_dir)
    times: list[float] = []
    fronts: list[float] = []
    for csv_path in frames:
        time = frame_time(csv_path)
        if not math.isfinite(time):
            continue
        rows = physical_rows(load_rows(csv_path))
        frame_profile = pick_row_by_y(rows)
        if not frame_profile:
            continue
        front = front_position_from_w(frame_profile)
        if math.isfinite(front):
            times.append(time)
            fronts.append(front)

    if times and fronts:
        fig, axis = plt.subplots(figsize=(9, 5))
        axis.plot(times, fronts, marker="o", linewidth=1.8)
        axis.set_xlabel("t [μs]")
        axis.set_ylabel("x_front [cm]")
        axis.set_title("CJ front tracking from W ≤ 0.5")
        axis.grid(True, alpha=0.3)
        annotate_units(fig)
        fig.tight_layout(rect=(0.04, 0.05, 0.98, 0.96))
        fig.savefig(plot_dir / "cj_front_tracking.png", dpi=160)
        plt.close(fig)


def build_grid(rows: list[dict[str, str]], field: str) -> tuple[list[float], list[float], list[list[float]]]:
    x_levels = sorted({f(row["x"]) for row in rows})
    y_levels = sorted({f(row["y"]) for row in rows})
    x_index = {value: idx for idx, value in enumerate(x_levels)}
    y_index = {value: idx for idx, value in enumerate(y_levels)}
    grid = [[math.nan for _ in x_levels] for _ in y_levels]
    for row in rows:
        x = f(row["x"])
        y = f(row["y"])
        value = f(row[field])
        if math.isfinite(x) and math.isfinite(y) and math.isfinite(value):
            grid[y_index[y]][x_index[x]] = value
    return x_levels, y_levels, grid


def plot_corner(case_dir: Path) -> None:
    plot_dir = ensure_plot_dir(case_dir)
    final_rows = physical_rows(load_rows(case_dir / "final_results.csv"))
    if not final_rows:
        return

    x_levels, y_levels, p_grid = build_grid(final_rows, "p")
    _, _, w_grid = build_grid(final_rows, "w")

    extent = [min(x_levels), max(x_levels), min(y_levels), max(y_levels)]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharex=True, sharey=True)
    images = [
        (axes[0], p_grid, "Pressure P [Mbar]", "inferno"),
        (axes[1], w_grid, "Unreacted mass fraction W [-]", "viridis_r"),
    ]
    for axis, grid, title, cmap in images:
        image = axis.imshow(grid, extent=extent, origin="lower", aspect="auto", cmap=cmap)
        axis.set_title(title)
        axis.set_xlabel("x [cm]")
        axis.grid(False)
        fig.colorbar(image, ax=axis, shrink=0.9)
    axes[0].set_ylabel("y [cm]")
    fig.suptitle("Corner turning, final state maps", fontsize=14)
    annotate_units(fig)
    fig.tight_layout(rect=(0.03, 0.05, 0.98, 0.95))
    fig.savefig(plot_dir / "corner_maps.png", dpi=160)
    plt.close(fig)

    frames = sorted_frames(case_dir)
    if not frames:
        return

    arrival_rows = physical_rows(load_rows(frames[-1]))
    x_levels, y_levels, arrival_grid = build_grid(arrival_rows, "p")
    for y_idx, _ in enumerate(y_levels):
        for x_idx, _ in enumerate(x_levels):
            arrival_grid[y_idx][x_idx] = math.nan

    x_index = {value: idx for idx, value in enumerate(x_levels)}
    y_index = {value: idx for idx, value in enumerate(y_levels)}

    for csv_path in frames:
        time = frame_time(csv_path)
        if not math.isfinite(time):
            continue
        for row in physical_rows(load_rows(csv_path)):
            x = f(row["x"])
            y = f(row["y"])
            w = f(row["w"])
            x_idx = x_index.get(x)
            y_idx = y_index.get(y)
            if x_idx is None or y_idx is None:
                continue
            if math.isfinite(w) and w <= 0.5 and not math.isfinite(arrival_grid[y_idx][x_idx]):
                arrival_grid[y_idx][x_idx] = time

    fig, axis = plt.subplots(figsize=(7, 5))
    image = axis.imshow(arrival_grid, extent=extent, origin="lower", aspect="auto", cmap="magma_r")
    axis.set_title("Corner turning arrival time from W ≤ 0.5")
    axis.set_xlabel("x [cm]")
    axis.set_ylabel("y [cm]")
    fig.colorbar(image, ax=axis, shrink=0.9, label="Arrival time [μs]")
    annotate_units(fig)
    fig.tight_layout(rect=(0.04, 0.05, 0.98, 0.96))
    fig.savefig(plot_dir / "corner_arrival_time.png", dpi=160)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Create presentation-ready plots for Mader 2DE test outputs")
    parser.add_argument("--case", required=True, choices=["sod", "cj", "corner"])
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()

    case_dir = Path(args.output_root).resolve() / args.case
    if not case_dir.exists():
        raise SystemExit(f"Missing case output directory: {case_dir}")

    if args.case == "sod":
        plot_sod(case_dir)
    elif args.case == "cj":
        plot_cj(case_dir)
    else:
        plot_corner(case_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
