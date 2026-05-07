#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from corner_radius_study import f, frame_time, infer_spacing, is_solid, load_rows, radius_tag, sorted_frames


@dataclass(frozen=True)
class CaseData:
    label: str
    radius: float
    case_dir: Path
    final_rows: list[dict[str, str]]
    frames: list[Path]
    p0: float
    p_ref: float
    p_threshold: float
    dx: float
    dy: float
    analysis_keys: set[tuple[str, str]]
    max_pressure: dict[tuple[str, str], float]
    arrival_time: dict[tuple[str, str], float]


FRONT_BANDS = {
    "top": (3.15, 3.95),
    "middle": (1.60, 2.40),
    "bottom": (0.20, 1.00),
}


def parse_radii(raw: str | None) -> list[float] | None:
    if raw is None:
        return None
    values = []
    for item in raw.split(","):
        item = item.strip()
        if item:
            values.append(float(item))
    return values


def parse_thresholds(raw: str) -> list[float]:
    values = [float(item.strip()) for item in raw.split(",") if item.strip()]
    if not values:
        raise SystemExit("At least one W threshold is required.")
    return values


def radius_from_tag(tag: str) -> float | None:
    if not tag.startswith("R"):
        return None
    text = tag[1:].replace("m", "-").replace("p", ".")
    try:
        return float(text)
    except ValueError:
        return None


def read_summary_radii(output_root: Path) -> list[float]:
    path = output_root / "summary.csv"
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return [float(row["R_corner_cm"]) for row in csv.DictReader(handle)]


def available_cases(output_root: Path, selected_radii: list[float] | None) -> list[tuple[float, Path]]:
    radii = read_summary_radii(output_root)
    if not radii:
        for child in output_root.iterdir():
            if child.is_dir():
                radius = radius_from_tag(child.name)
                if radius is not None:
                    radii.append(radius)
    if selected_radii is not None:
        wanted = set(selected_radii)
        radii = [radius for radius in radii if radius in wanted]
    cases = []
    for radius in sorted(set(radii)):
        case_dir = output_root / radius_tag(radius)
        if case_dir.exists():
            cases.append((radius, case_dir))
    return cases


def infer_grid_label(output_root: Path, explicit: str | None) -> str:
    if explicit:
        return explicit
    match = re.search(r"(\d+)x(\d+)", output_root.name)
    if match:
        return match.group(0)
    first_summary = read_summary_radii(output_root)
    if first_summary:
        case_dir = output_root / radius_tag(first_summary[0])
        final = case_dir / "final_results.csv"
        if final.exists():
            rows = load_rows(final)
            xs = {row["x"] for row in rows}
            ys = {row["y"] for row in rows}
            return f"{len(xs)}x{len(ys)}"
    return output_root.name


def load_case(
    label: str,
    radius: float,
    case_dir: Path,
    x_analysis_min: float,
    y_analysis_max: float,
    arrival_fraction: float,
) -> CaseData:
    initial = case_dir / "step_0_initial.csv"
    final = case_dir / "final_results.csv"
    if not initial.exists() or not final.exists():
        raise SystemExit(f"Missing initial/final CSV in {case_dir}")

    initial_rows = load_rows(initial)
    p_values = [f(row["p"]) for row in initial_rows if not is_solid(row) and math.isfinite(f(row["p"]))]
    if not p_values:
        raise SystemExit(f"Could not infer P_ref from {initial}")
    p0 = min(p_values)
    p_ref = max(p_values)
    p_threshold = p0 + arrival_fraction * (p_ref - p0)

    frames = sorted_frames(case_dir)
    max_pressure: dict[tuple[str, str], float] = {}
    arrival_time: dict[tuple[str, str], float] = {}
    for csv_path in frames:
        time = frame_time(csv_path)
        for row in load_rows(csv_path):
            if is_solid(row):
                continue
            pressure = f(row["p"])
            if not math.isfinite(pressure):
                continue
            key = (row["x"], row["y"])
            max_pressure[key] = max(max_pressure.get(key, -math.inf), pressure)
            if math.isfinite(time) and pressure >= p_threshold and key not in arrival_time:
                arrival_time[key] = time

    final_rows = load_rows(final)
    dx, dy = infer_spacing(final_rows)
    analysis_keys: set[tuple[str, str]] = set()
    for row in final_rows:
        if is_solid(row):
            continue
        x = f(row["x"])
        y = f(row["y"])
        if not math.isfinite(x) or not math.isfinite(y):
            continue
        if x + 1e-12 < x_analysis_min or y - 1e-12 > y_analysis_max:
            continue
        key = (row["x"], row["y"])
        if max_pressure.get(key, -math.inf) >= p_threshold:
            analysis_keys.add(key)

    return CaseData(label, radius, case_dir, final_rows, frames, p0, p_ref, p_threshold, dx, dy, analysis_keys, max_pressure, arrival_time)


def write_threshold_sensitivity(out_dir: Path, cases: list[CaseData], thresholds: list[float], epsilon: float) -> Path:
    rows = []
    for case in cases:
        cell_area = case.dx * case.dy
        s_analysis = len(case.analysis_keys) * cell_area
        box_cells = sum(
            1 for row in case.final_rows
            if not is_solid(row)
            and math.isfinite(f(row["x"]))
            and math.isfinite(f(row["y"]))
            and f(row["x"]) >= 2.0
            and f(row["y"]) <= 3.0
        )
        s_box = box_cells * cell_area
        arrived_fraction = len(case.analysis_keys) / box_cells if box_cells > 0 else math.nan
        final_by_key = {(row["x"], row["y"]): row for row in case.final_rows}
        for threshold in thresholds:
            dead_cells = 0
            for key in case.analysis_keys:
                row = final_by_key.get(key)
                if row is not None and f(row["w"]) > threshold:
                    dead_cells += 1
            s_dead = dead_cells * cell_area
            rows.append({
                "grid": case.label,
                "R_corner_cm": case.radius,
                "W_thr": threshold,
                "S_dead_cm2": s_dead,
                "S_hat": s_dead / s_analysis if s_analysis > 0.0 else math.nan,
                "S_analysis_cm2": s_analysis,
                "S_box_cm2": s_box,
                "arrived_fraction": arrived_fraction,
                "p0": case.p0,
                "P_ref": case.p_ref,
                "p_arrival_threshold": case.p_threshold,
                "dead_cells": dead_cells,
                "arrived_cells": len(case.analysis_keys),
            })

    path = out_dir / "threshold_sensitivity.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    critical_rows = []
    for label in sorted({row["grid"] for row in rows}):
        for threshold in thresholds:
            data = sorted(
                (row for row in rows if row["grid"] == label and abs(row["W_thr"] - threshold) < 1e-12),
                key=lambda item: item["R_corner_cm"],
            )
            critical = next((row["R_corner_cm"] for row in data if row["S_hat"] < epsilon), math.nan)
            text = f"{critical:g}" if math.isfinite(critical) else f">{max(row['R_corner_cm'] for row in data):g}"
            critical_rows.append({
                "grid": label,
                "W_thr": threshold,
                "epsilon": epsilon,
                "R_cr_cm": text,
            })

    with (out_dir / "critical_by_threshold.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["grid", "W_thr", "epsilon", "R_cr_cm"])
        writer.writeheader()
        writer.writerows(critical_rows)

    fig, axis = plt.subplots(figsize=(8.5, 5.0))
    for threshold in thresholds:
        subset = [row for row in rows if abs(row["W_thr"] - threshold) < 1e-12]
        labels = sorted({row["grid"] for row in subset})
        for label in labels:
            data = sorted((row for row in subset if row["grid"] == label), key=lambda item: item["R_corner_cm"])
            axis.plot([row["R_corner_cm"] for row in data],
                      [row["S_hat"] for row in data],
                      marker="o",
                      linewidth=1.5,
                      label=f"{label}, W>{threshold:g}")
    axis.set_xlabel("R_corner [cm]")
    axis.set_ylabel("S_dead / S_analysis [-]")
    axis.set_title("Dead-zone sensitivity to W threshold")
    axis.grid(True, alpha=0.3)
    axis.legend(fontsize=8, ncols=2)
    fig.tight_layout()
    fig.savefig(out_dir / "threshold_sensitivity.png", dpi=170)
    plt.close(fig)

    finest_label = max(sorted({row["grid"] for row in rows}), key=lambda label: int(label.split("x", 1)[0]) if "x" in label and label.split("x", 1)[0].isdigit() else 0)
    fig, axis = plt.subplots(figsize=(7.6, 4.8))
    for threshold in thresholds:
        data = sorted(
            (row for row in rows if row["grid"] == finest_label and abs(row["W_thr"] - threshold) < 1e-12),
            key=lambda item: item["R_corner_cm"],
        )
        axis.plot([row["R_corner_cm"] for row in data],
                  [row["S_hat"] for row in data],
                  marker="o",
                  linewidth=1.8,
                  label=f"W_thr={threshold:g}")
    axis.axhline(epsilon, color="black", linestyle="--", linewidth=1.2, label=f"epsilon={epsilon:g}")
    axis.set_xlabel("R_corner [cm]")
    axis.set_ylabel("S_dead / S_analysis [-]")
    axis.set_title(f"Threshold sensitivity on {finest_label}")
    axis.grid(True, alpha=0.3)
    axis.legend(fontsize=9, ncols=2)
    fig.tight_layout()
    fig.savefig(out_dir / "threshold_sensitivity_fine.png", dpi=170)
    plt.close(fig)
    return path


def front_position(rows: list[dict[str, str]], p_threshold: float, band: tuple[float, float]) -> float:
    y_min, y_max = band
    x_front = math.nan
    for row in rows:
        if is_solid(row):
            continue
        x = f(row["x"])
        y = f(row["y"])
        p = f(row["p"])
        if not math.isfinite(x) or not math.isfinite(y) or not math.isfinite(p):
            continue
        if x < 2.0 or y < y_min or y > y_max:
            continue
        if p >= p_threshold:
            x_front = x if not math.isfinite(x_front) else max(x_front, x)
    return x_front


def linear_fit(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    if len(xs) < 2:
        return math.nan, math.nan, math.nan
    x_mean = sum(xs) / len(xs)
    y_mean = sum(ys) / len(ys)
    sxx = sum((x - x_mean) ** 2 for x in xs)
    if sxx <= 0.0:
        return math.nan, math.nan, math.nan
    sxy = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    slope = sxy / sxx
    intercept = y_mean - slope * x_mean
    ss_tot = sum((y - y_mean) ** 2 for y in ys)
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, ys))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else math.nan
    return slope, intercept, r2


def write_front_metrics(out_dir: Path, cases: list[CaseData], arrival_fraction: float) -> tuple[Path, Path]:
    trajectory_rows = []
    metric_rows = []
    for case in cases:
        p_threshold = case.p_threshold
        by_band: dict[str, list[tuple[float, float]]] = {name: [] for name in FRONT_BANDS}
        for csv_path in case.frames:
            time = frame_time(csv_path)
            if not math.isfinite(time) or time < 1.5:
                continue
            rows = load_rows(csv_path)
            for name, band in FRONT_BANDS.items():
                x_front = front_position(rows, p_threshold, band)
                if math.isfinite(x_front):
                    by_band[name].append((time, x_front))
                    trajectory_rows.append({
                        "grid": case.label,
                        "R_corner_cm": case.radius,
                        "band": name,
                        "time_us": time,
                        "x_front_cm": x_front,
                    })

        final_positions = {
            name: values[-1][1] for name, values in by_band.items() if values
        }
        top = final_positions.get("top", math.nan)
        bottom = final_positions.get("bottom", math.nan)
        for name, values in by_band.items():
            fit_values = [(t, x) for t, x in values if 2.2 <= x <= 6.6]
            slope, intercept, r2 = linear_fit([t for t, _ in fit_values], [x for _, x in fit_values])
            metric_rows.append({
                "grid": case.label,
                "R_corner_cm": case.radius,
                "band": name,
                "front_speed_cm_per_us": slope,
                "fit_r2": r2,
                "last_x_front_cm": values[-1][1] if values else math.nan,
                "top_minus_bottom_lag_cm": top - bottom if name == "bottom" and math.isfinite(top) and math.isfinite(bottom) else math.nan,
                "points_used": len(fit_values),
            })

    trajectory_path = out_dir / "front_trajectories.csv"
    with trajectory_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["grid", "R_corner_cm", "band", "time_us", "x_front_cm"])
        writer.writeheader()
        writer.writerows(trajectory_rows)

    metrics_path = out_dir / "front_metrics.csv"
    with metrics_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metric_rows[0].keys()))
        writer.writeheader()
        writer.writerows(metric_rows)

    labels = sorted({row["grid"] for row in trajectory_rows})
    for label in labels:
        label_rows = [row for row in trajectory_rows if row["grid"] == label]
        radii = sorted({row["R_corner_cm"] for row in label_rows})
        fig, axes = plt.subplots(len(radii), 1, figsize=(8.0, max(3.0, 2.2 * len(radii))), sharex=True, sharey=True)
        if len(radii) == 1:
            axes = [axes]
        for axis, radius in zip(axes, radii):
            for band in FRONT_BANDS:
                data = [row for row in label_rows if row["R_corner_cm"] == radius and row["band"] == band]
                axis.plot([row["time_us"] for row in data], [row["x_front_cm"] for row in data], marker=".", linewidth=1.2, label=band)
            axis.set_title(f"{label}, R = {radius:g} cm")
            axis.set_ylabel("x_front [cm]")
            axis.grid(True, alpha=0.3)
        axes[-1].set_xlabel("t [us]")
        axes[0].legend(fontsize=8, ncols=3)
        fig.tight_layout()
        fig.savefig(out_dir / f"front_trajectories_{label}.png", dpi=170)
        plt.close(fig)

    return trajectory_path, metrics_path


def write_arrival_similarity(out_dir: Path, cases: list[CaseData]) -> Path:
    rows = []
    for case in cases:
        times = []
        radii = []
        for row in case.final_rows:
            if is_solid(row):
                continue
            x = f(row["x"])
            y = f(row["y"])
            if not math.isfinite(x) or not math.isfinite(y):
                continue
            key = (row["x"], row["y"])
            t = case.arrival_time.get(key, math.nan)
            if not math.isfinite(t) or x < 2.0 or y > 3.0:
                continue
            r = math.hypot(x - 2.0, y - 3.0)
            if r < 0.25:
                continue
            times.append(t)
            radii.append(r)
        speed, intercept, r2 = linear_fit(times, radii)
        rows.append({
            "grid": case.label,
            "R_corner_cm": case.radius,
            "cylindrical_proxy_speed_cm_per_us": speed,
            "intercept_cm": intercept,
            "arrival_r2": r2,
            "points_used": len(times),
        })

    path = out_dir / "arrival_similarity.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fig, axis = plt.subplots(figsize=(7.5, 4.8))
    for label in sorted({row["grid"] for row in rows}):
        data = sorted((row for row in rows if row["grid"] == label), key=lambda item: item["R_corner_cm"])
        axis.plot([row["R_corner_cm"] for row in data], [row["arrival_r2"] for row in data], marker="o", label=label)
    axis.set_xlabel("R_corner [cm]")
    axis.set_ylabel("R^2 of r(corner) vs arrival time")
    axis.set_title("Cylindrical self-similarity proxy for diffracted front")
    axis.set_ylim(0.0, 1.02)
    axis.grid(True, alpha=0.3)
    axis.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "arrival_similarity.png", dpi=170)
    plt.close(fig)
    return path


def write_grid_convergence(out_dir: Path, cases: list[CaseData], w_threshold: float) -> Path:
    rows = []
    for case in cases:
        h = math.sqrt(case.dx * case.dy)
        s_analysis = len(case.analysis_keys) * case.dx * case.dy
        final_by_key = {(row["x"], row["y"]): row for row in case.final_rows}
        dead_cells = sum(
            1 for key in case.analysis_keys
            if key in final_by_key and f(final_by_key[key]["w"]) > w_threshold
        )
        s_dead = dead_cells * case.dx * case.dy
        rows.append({
            "grid": case.label,
            "R_corner_cm": case.radius,
            "h_cm": h,
            "S_dead_cm2": s_dead,
            "S_hat": s_dead / s_analysis if s_analysis > 0.0 else math.nan,
        })

    path = out_dir / "grid_convergence_extended.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fig, axis = plt.subplots(figsize=(7.5, 4.8))
    for radius in sorted({row["R_corner_cm"] for row in rows}):
        data = sorted((row for row in rows if row["R_corner_cm"] == radius), key=lambda item: item["h_cm"], reverse=True)
        axis.plot([row["h_cm"] for row in data], [row["S_hat"] for row in data], marker="o", label=f"R={radius:g} cm")
    axis.set_xlabel("effective grid spacing h [cm]")
    axis.set_ylabel("S_dead / S_analysis [-]")
    axis.set_title(f"Grid convergence of dead-zone metric, W>{w_threshold:g}")
    axis.grid(True, alpha=0.3)
    axis.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "grid_convergence_extended.png", dpi=170)
    plt.close(fig)
    return path


def grid_width(label: str) -> int:
    if "x" not in label:
        return 0
    head = label.split("x", 1)[0]
    return int(head) if head.isdigit() else 0


def write_radius_curves(out_dir: Path, cases: list[CaseData], w_threshold: float, epsilon: float) -> Path:
    rows = []
    for case in cases:
        s_analysis = len(case.analysis_keys) * case.dx * case.dy
        box_cells = sum(
            1 for row in case.final_rows
            if not is_solid(row)
            and math.isfinite(f(row["x"]))
            and math.isfinite(f(row["y"]))
            and f(row["x"]) >= 2.0
            and f(row["y"]) <= 3.0
        )
        s_box = box_cells * case.dx * case.dy
        final_by_key = {(row["x"], row["y"]): row for row in case.final_rows}
        dead_cells = sum(
            1 for key in case.analysis_keys
            if key in final_by_key and f(final_by_key[key]["w"]) > w_threshold
        )
        s_dead = dead_cells * case.dx * case.dy
        rows.append({
            "grid": case.label,
            "R_corner_cm": case.radius,
            "S_dead_cm2": s_dead,
            "S_hat": s_dead / s_analysis if s_analysis > 0.0 else math.nan,
            "S_analysis_cm2": s_analysis,
            "S_box_cm2": s_box,
            "arrived_fraction": len(case.analysis_keys) / box_cells if box_cells > 0 else math.nan,
        })

    path = out_dir / "s_hat_by_radius.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    labels = sorted({row["grid"] for row in rows}, key=grid_width)
    for filename, min_width in (("s_hat_by_radius.png", 0), ("s_hat_by_radius_fine.png", 560)):
        fig, axis = plt.subplots(figsize=(7.6, 4.8))
        for label in labels:
            if grid_width(label) < min_width:
                continue
            data = sorted((row for row in rows if row["grid"] == label), key=lambda item: item["R_corner_cm"])
            axis.plot(
                [row["R_corner_cm"] for row in data],
                [row["S_hat"] for row in data],
                marker="o",
                linewidth=1.8,
                label=label,
            )
        axis.axhline(epsilon, color="black", linestyle="--", linewidth=1.2, label=f"epsilon={epsilon:g}")
        axis.set_xlabel("R_corner [cm]")
        axis.set_ylabel("S_dead / S_analysis [-]")
        title = "Relative dead-zone area vs radius"
        if min_width:
            title += ", fine grids"
        axis.set_title(title)
        axis.grid(True, alpha=0.3)
        axis.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(out_dir / filename, dpi=170)
        plt.close(fig)

    for filename, min_width in (("arrived_fraction_by_radius.png", 0), ("arrived_fraction_by_radius_fine.png", 560)):
        fig, axis = plt.subplots(figsize=(7.6, 4.8))
        for label in labels:
            if grid_width(label) < min_width:
                continue
            data = sorted((row for row in rows if row["grid"] == label), key=lambda item: item["R_corner_cm"])
            axis.plot(
                [row["R_corner_cm"] for row in data],
                [row["arrived_fraction"] for row in data],
                marker="o",
                linewidth=1.8,
                label=label,
            )
        axis.set_xlabel("R_corner [cm]")
        axis.set_ylabel("S_analysis / S_box [-]")
        title = "Arrived fraction of downstream box"
        if min_width:
            title += ", fine grids"
        axis.set_title(title)
        axis.grid(True, alpha=0.3)
        axis.legend(fontsize=9)
        fig.tight_layout()
        fig.savefig(out_dir / filename, dpi=170)
        plt.close(fig)
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description="Post-process rounded-corner detonation validation metrics.")
    parser.add_argument("--output-root", action="append", required=True, help="Study output root. May be repeated.")
    parser.add_argument("--label", action="append", default=None, help="Optional label for each output root.")
    parser.add_argument("--radii", default=None, help="Optional comma-separated radii filter.")
    parser.add_argument("--thresholds", default="0.02,0.05,0.10,0.20")
    parser.add_argument("--out-dir", default="tests/mader_2de/report/validation")
    parser.add_argument("--x-analysis-min", type=float, default=2.0)
    parser.add_argument("--y-analysis-max", type=float, default=3.0)
    parser.add_argument("--arrival-fraction", type=float, default=0.01)
    parser.add_argument("--grid-w-threshold", type=float, default=0.05)
    parser.add_argument("--epsilon", type=float, default=0.01)
    args = parser.parse_args()

    roots = [Path(item).resolve() for item in args.output_root]
    labels = args.label or []
    if labels and len(labels) != len(roots):
        raise SystemExit("--label must be supplied once per --output-root")

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    selected = parse_radii(args.radii)
    thresholds = parse_thresholds(args.thresholds)

    cases: list[CaseData] = []
    for idx, root in enumerate(roots):
        label = labels[idx] if labels else infer_grid_label(root, None)
        for radius, case_dir in available_cases(root, selected):
            print(f"[load] {label} R={radius:g}: {case_dir}", flush=True)
            cases.append(load_case(label, radius, case_dir, args.x_analysis_min, args.y_analysis_max, args.arrival_fraction))

    if not cases:
        raise SystemExit("No cases found.")

    print(f"[write] {write_threshold_sensitivity(out_dir, cases, thresholds, args.epsilon)}")
    trajectory_path, metrics_path = write_front_metrics(out_dir, cases, args.arrival_fraction)
    print(f"[write] {trajectory_path}")
    print(f"[write] {metrics_path}")
    print(f"[write] {write_arrival_similarity(out_dir, cases)}")
    print(f"[write] {write_grid_convergence(out_dir, cases, args.grid_w_threshold)}")
    print(f"[write] {write_radius_curves(out_dir, cases, args.grid_w_threshold, args.epsilon)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
