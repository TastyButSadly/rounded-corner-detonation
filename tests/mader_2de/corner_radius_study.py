#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_RADII = "0,0.25,0.5,1.0,1.5,2.0"


@dataclass(frozen=True)
class StudyResult:
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


def resolve_path(root: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = root / path
    return path.resolve()


def resolve_solver(root: Path, override: str | None) -> Path:
    if override:
        solver = resolve_path(root, override)
        if not solver.exists():
            raise SystemExit(f"Solver not found: {solver}")
        return solver

    for candidate in (
        root / "build" / "GodunovSolver",
        root / "build-codex" / "GodunovSolver",
        root / "build-mpi" / "GodunovSolver",
    ):
        if candidate.exists():
            return candidate.resolve()
    raise SystemExit("Solver not found. Use --solver to specify the executable.")


def parse_radii(raw: str) -> list[float]:
    radii = []
    for item in raw.split(","):
        item = item.strip()
        if not item:
            continue
        radius = float(item)
        if radius < 0.0:
            raise SystemExit("Radii must be non-negative.")
        radii.append(radius)
    if not radii:
        raise SystemExit("At least one radius is required.")
    return radii


def radius_tag(radius: float) -> str:
    text = f"{radius:g}".replace("-", "m").replace(".", "p")
    return f"R{text}"


def format_number(value: float) -> str:
    return f"{value:.12g}"


def key_from_line(line: str) -> str | None:
    body = line.split(";", 1)[0].split("#", 1)[0].strip()
    if "=" not in body:
        return None
    return body.split("=", 1)[0].strip()


def update_config_text(
    text: str,
    system_values: dict[str, str],
    block_values: dict[str, dict[str, str]],
) -> str:
    lines = text.splitlines()
    out: list[str] = []
    current_section = ""
    current_block = ""
    pending_system = dict(system_values)

    def flush_system_values() -> None:
        nonlocal pending_system
        for key, value in pending_system.items():
            out.append(f"{key} = {value}")
        pending_system = {}

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            if current_section == "SYSTEM":
                flush_system_values()
            current_section = stripped[1:-1].strip()
            current_block = ""
            out.append(line)
            continue

        block_match = re.fullmatch(r"BLOCK_\d+", stripped)
        if current_section == "INITIAL_CONDITIONS" and block_match:
            current_block = stripped
            out.append(line)
            continue

        key = key_from_line(line)
        if current_section == "SYSTEM" and key in system_values:
            out.append(f"{key} = {system_values[key]}")
            pending_system.pop(key, None)
            continue

        if current_section == "INITIAL_CONDITIONS" and current_block in block_values:
            updates = block_values[current_block]
            if key in updates:
                out.append(f"{key} = {updates[key]}")
                continue

        out.append(line)

    if current_section == "SYSTEM":
        flush_system_values()
    return "\n".join(out) + "\n"


def generate_config(
    base_config: Path,
    output_config: Path,
    radius: float,
    tag: str,
    nx: int | None,
    ny: int | None,
    tmax: float | None,
    dt_out: float | None,
) -> None:
    system_values = {
        "case_name": f"mader_2de_corner_radius_{tag}",
        "corner_radius": format_number(radius),
    }
    block_values: dict[str, dict[str, str]] = {}

    if nx is not None:
        system_values["Nx"] = str(nx)
        block_values.setdefault("BLOCK_2", {})["x_end"] = str(nx - 1)
    if ny is not None:
        system_values["Ny"] = str(ny)
        block_values.setdefault("BLOCK_1", {})["y_end"] = str(ny - 1)
        block_values.setdefault("BLOCK_2", {})["y_end"] = str(ny - 1)
    if tmax is not None:
        system_values["tmax"] = format_number(tmax)
    if dt_out is not None:
        system_values["dt_out"] = format_number(dt_out)

    config_text = update_config_text(base_config.read_text(encoding="utf-8"), system_values, block_values)
    output_config.parent.mkdir(parents=True, exist_ok=True)
    output_config.write_text(config_text, encoding="utf-8")


def run_solver(solver: Path, config: Path, case_dir: Path, root: Path, keep_existing: bool) -> None:
    if case_dir.exists() and not keep_existing:
        shutil.rmtree(case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)

    print(f"[run] {case_dir.name}: {solver} {config} {case_dir}", flush=True)
    proc = subprocess.run([str(solver), str(config), str(case_dir)], cwd=str(root), env=os.environ.copy())
    if proc.returncode != 0:
        raise SystemExit(f"Solver failed for {case_dir.name} with code {proc.returncode}")


def f(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return math.nan


def is_solid(row: dict[str, str]) -> bool:
    return int(float(row["is_solid"])) != 0


def load_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def frame_time(csv_path: Path) -> float:
    if csv_path.name == "step_0_initial.csv":
        return 0.0
    if csv_path.name == "final_results.csv":
        return math.inf
    marker = "_time_"
    if marker in csv_path.name:
        try:
            return float(csv_path.name.split(marker, 1)[1].rsplit(".csv", 1)[0])
        except ValueError:
            return math.nan
    return math.nan


def sorted_frames(case_dir: Path) -> list[Path]:
    frames: list[Path] = []
    initial = case_dir / "step_0_initial.csv"
    final = case_dir / "final_results.csv"
    if initial.exists():
        frames.append(initial)
    frames.extend(sorted(case_dir.glob("step_*_time_*.csv"), key=frame_time))
    if final.exists():
        frames.append(final)
    return frames


def infer_spacing(rows: list[dict[str, str]]) -> tuple[float, float]:
    x_values = sorted({f(row["x"]) for row in rows if math.isfinite(f(row["x"]))})
    y_values = sorted({f(row["y"]) for row in rows if math.isfinite(f(row["y"]))})
    dx_values = [b - a for a, b in zip(x_values, x_values[1:]) if b > a]
    dy_values = [b - a for a, b in zip(y_values, y_values[1:]) if b > a]
    if not dx_values or not dy_values:
        raise SystemExit("Could not infer grid spacing from final_results.csv")
    return min(dx_values), min(dy_values)


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


def plot_field(rows: list[dict[str, str]], field: str, title: str, path: Path, cmap: str, vmin: float | None, vmax: float | None) -> None:
    x_levels, y_levels, grid = build_grid(rows, field)
    extent = [min(x_levels), max(x_levels), min(y_levels), max(y_levels)]
    fig, axis = plt.subplots(figsize=(8.5, 4.8))
    image = axis.imshow(grid, extent=extent, origin="lower", aspect="equal", cmap=cmap, vmin=vmin, vmax=vmax)
    axis.set_title(title)
    axis.set_xlabel("x [cm]")
    axis.set_ylabel("y [cm]")
    fig.colorbar(image, ax=axis, shrink=0.9)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def plot_dead_map(
    rows: list[dict[str, str]],
    analysis_keys: set[tuple[str, str]],
    dead_keys: set[tuple[str, str]],
    title: str,
    path: Path,
) -> None:
    x_levels = sorted({f(row["x"]) for row in rows})
    y_levels = sorted({f(row["y"]) for row in rows})
    x_index = {value: idx for idx, value in enumerate(x_levels)}
    y_index = {value: idx for idx, value in enumerate(y_levels)}
    grid = [[math.nan for _ in x_levels] for _ in y_levels]

    for row in rows:
        x = f(row["x"])
        y = f(row["y"])
        if not math.isfinite(x) or not math.isfinite(y) or is_solid(row):
            continue
        key = (row["x"], row["y"])
        if key in analysis_keys:
            grid[y_index[y]][x_index[x]] = 1.0 if key in dead_keys else 0.0

    cmap = plt.get_cmap("gray_r").copy()
    cmap.set_bad("#dddddd")
    extent = [min(x_levels), max(x_levels), min(y_levels), max(y_levels)]
    fig, axis = plt.subplots(figsize=(8.5, 4.8))
    image = axis.imshow(grid, extent=extent, origin="lower", aspect="equal", cmap=cmap, vmin=0.0, vmax=1.0)
    axis.set_title(title)
    axis.set_xlabel("x [cm]")
    axis.set_ylabel("y [cm]")
    fig.colorbar(image, ax=axis, shrink=0.9, label="dead = 1")
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def analyze_case(
    case_dir: Path,
    radius: float,
    x_analysis_min: float,
    y_analysis_max: float,
    w_threshold: float,
) -> StudyResult:
    frames = sorted_frames(case_dir)
    if not frames:
        raise SystemExit(f"No CSV frames found in {case_dir}")

    initial = case_dir / "step_0_initial.csv"
    final = case_dir / "final_results.csv"
    if not initial.exists() or not final.exists():
        raise SystemExit(f"Missing initial or final CSV in {case_dir}")

    initial_rows = load_rows(initial)
    initial_pressures = [f(row["p"]) for row in initial_rows if not is_solid(row) and math.isfinite(f(row["p"]))]
    if not initial_pressures:
        raise SystemExit(f"Could not determine P_ref in {initial}")
    p0 = min(initial_pressures)
    p_ref = max(initial_pressures)
    p_arrival_threshold = p0 + 0.01 * (p_ref - p0)

    max_pressure: dict[tuple[str, str], float] = {}
    for csv_path in frames:
        for row in load_rows(csv_path):
            if is_solid(row):
                continue
            pressure = f(row["p"])
            if not math.isfinite(pressure):
                continue
            key = (row["x"], row["y"])
            max_pressure[key] = max(max_pressure.get(key, -math.inf), pressure)

    final_rows = load_rows(final)
    dx, dy = infer_spacing(final_rows)
    cell_area = dx * dy
    analysis_keys: set[tuple[str, str]] = set()
    dead_keys: set[tuple[str, str]] = set()
    box_cells = 0

    for row in final_rows:
        if is_solid(row):
            continue
        x = f(row["x"])
        y = f(row["y"])
        if not math.isfinite(x) or not math.isfinite(y):
            continue
        if x + 1e-12 < x_analysis_min or y - 1e-12 > y_analysis_max:
            continue

        box_cells += 1
        key = (row["x"], row["y"])
        if max_pressure.get(key, -math.inf) < p_arrival_threshold:
            continue

        analysis_keys.add(key)
        w_value = f(row["w"])
        if math.isfinite(w_value) and w_value > w_threshold:
            dead_keys.add(key)

    arrived_cells = len(analysis_keys)
    dead_cells = len(dead_keys)
    s_analysis = arrived_cells * cell_area
    s_box = box_cells * cell_area
    s_dead = dead_cells * cell_area
    s_hat = s_dead / s_analysis if s_analysis > 0.0 else math.nan
    arrived_fraction = arrived_cells / box_cells if box_cells > 0 else math.nan

    plot_dir = case_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)
    tag = radius_tag(radius)
    plot_field(final_rows, "p", f"Final pressure, R = {radius:g} cm", plot_dir / f"{tag}_final_pressure.png", "inferno", None, None)
    plot_field(final_rows, "w", f"Final W, R = {radius:g} cm", plot_dir / f"{tag}_final_w.png", "viridis_r", 0.0, 1.0)
    plot_dead_map(final_rows, analysis_keys, dead_keys, f"Dead-zone map, R = {radius:g} cm", plot_dir / f"{tag}_dead_zone.png")

    return StudyResult(radius, s_dead, s_hat, s_analysis, s_box, arrived_fraction, arrived_cells, dead_cells)


def write_summary(output_root: Path, results: list[StudyResult], epsilon: float) -> None:
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
        max_radius = max(result.radius for result in results)
        text = f"R_cr > {max_radius:g} cm for epsilon = {epsilon:g}\n"
    else:
        text = f"R_cr = {critical:g} cm for epsilon = {epsilon:g}\n"
    (output_root / "critical_radius.txt").write_text(text, encoding="utf-8")


def result_is_critical(result: StudyResult, epsilon: float) -> bool:
    return math.isfinite(result.s_hat) and result.s_hat < epsilon


def plot_summary(output_root: Path, results: list[StudyResult]) -> None:
    sorted_results = sorted(results, key=lambda item: item.radius)
    radii = [result.radius for result in sorted_results]
    s_dead = [result.s_dead for result in sorted_results]
    s_hat = [result.s_hat for result in sorted_results]

    plot_dir = output_root / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    fig, axis = plt.subplots(figsize=(7.5, 4.8))
    axis.plot(radii, s_dead, marker="o", linewidth=1.8)
    axis.set_xlabel("R_corner [cm]")
    axis.set_ylabel("S_dead [cm^2]")
    axis.set_title("Dead-zone area vs rounded-corner radius")
    axis.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(plot_dir / "S_dead_vs_R.png", dpi=160)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7.5, 4.8))
    axis.plot(radii, s_hat, marker="o", linewidth=1.8)
    axis.set_xlabel("R_corner [cm]")
    axis.set_ylabel("S_dead / S_analysis [-]")
    axis.set_title("Relative dead-zone area vs rounded-corner radius")
    axis.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(plot_dir / "S_hat_vs_R.png", dpi=160)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run and analyze the Mader 2DE rounded-corner detonation study.")
    parser.add_argument("--solver", default=None, help="Path to GodunovSolver. Defaults to build, build-codex, or build-mpi.")
    parser.add_argument("--base-config", default="tests/mader_2de/configs/corner_radius_study_base.ini")
    parser.add_argument("--output-root", default="tests/mader_2de/output_corner_radius_study")
    parser.add_argument("--radii", default=DEFAULT_RADII)
    parser.add_argument("--nx", type=int, default=None, help="Optional Nx override for generated configs.")
    parser.add_argument("--ny", type=int, default=None, help="Optional Ny override for generated configs.")
    parser.add_argument("--tmax", type=float, default=None, help="Optional tmax override for generated configs.")
    parser.add_argument("--dt-out", type=float, default=None, help="Optional dt_out override for generated configs.")
    parser.add_argument("--x-analysis-min", type=float, default=2.0)
    parser.add_argument("--y-analysis-max", type=float, default=3.0)
    parser.add_argument("--w-threshold", type=float, default=0.05)
    parser.add_argument("--epsilon", type=float, default=0.01)
    parser.add_argument("--skip-run", action="store_true", help="Analyze existing per-radius outputs without running the solver.")
    parser.add_argument("--keep-existing", action="store_true", help="Do not delete an existing case directory before a solver run.")
    args = parser.parse_args()

    root = repo_root()
    base_config = resolve_path(root, args.base_config)
    output_root = resolve_path(root, args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    if not base_config.exists():
        raise SystemExit(f"Base config not found: {base_config}")

    if (args.nx is None) != (args.ny is None):
        raise SystemExit("--nx and --ny must be provided together.")
    if args.nx is not None and (args.nx < 2 or args.ny < 2):
        raise SystemExit("--nx and --ny must be at least 2.")

    solver = None if args.skip_run else resolve_solver(root, args.solver)
    radii = parse_radii(args.radii)
    config_dir = output_root / "configs"
    results: list[StudyResult] = []

    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str((root / "tests" / "mader_2de" / ".mplcache").resolve()))

    for radius in radii:
        tag = radius_tag(radius)
        config_path = config_dir / f"corner_{tag}.ini"
        case_dir = output_root / tag
        generate_config(base_config, config_path, radius, tag, args.nx, args.ny, args.tmax, args.dt_out)
        if solver is not None:
            run_solver(solver, config_path, case_dir, root, args.keep_existing)
        result = analyze_case(case_dir, radius, args.x_analysis_min, args.y_analysis_max, args.w_threshold)
        results.append(result)
        print(
            f"[analyze] R={radius:g}: S_dead={result.s_dead:.6g}, "
            f"S_hat={result.s_hat:.6g}, dead_cells={result.dead_cells}, arrived_cells={result.arrived_cells}",
            flush=True,
        )

    write_summary(output_root, results, args.epsilon)
    plot_summary(output_root, results)
    print(f"[summary] wrote {output_root / 'summary.csv'}")
    print(f"[summary] wrote {output_root / 'critical_radius.txt'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
