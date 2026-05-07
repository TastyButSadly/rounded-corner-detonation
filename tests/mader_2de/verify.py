#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import statistics
from dataclasses import dataclass
from pathlib import Path


REQUIRED_COLUMNS = [
    "x",
    "y",
    "rho",
    "u",
    "v",
    "p",
    "w",
    "temperature",
    "is_solid",
]

SOD_EXPECTED = {
    "rho_left": 0.4263,
    "p_left": 0.3031,
    "u_left": 0.9274,
    "rho_right": 0.2656,
    "rho_far_right": 0.125,
    "p_right": 0.3031,
    "u_right": 0.9274,
    "contact_x": 0.685,
    "shock_x": 0.85,
}

CJ_EXPECTED = {
    "gamma": 3.0,
    "rho0": 1.84,
    "detonator_rho": 1.84 * 4.0 / 3.0,
    "detonator_cells": 5,
    "d_cj": 0.88,
}
CJ_EXPECTED["p_cj"] = CJ_EXPECTED["d_cj"] ** 2 * CJ_EXPECTED["rho0"] / (CJ_EXPECTED["gamma"] + 1.0)

CORNER_EXPECTED = {
    "gamma": 3.0,
    "rho0": 1.84,
    "d_cj": 0.88,
    "channel_exit_x": 2.0,
    "channel_top_y": 3.0,
}
CORNER_EXPECTED["p_cj"] = CORNER_EXPECTED["d_cj"] ** 2 * CORNER_EXPECTED["rho0"] / (CORNER_EXPECTED["gamma"] + 1.0)


@dataclass(frozen=True)
class CaseResult:
    name: str
    passed: bool
    message: str


def load_rows(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open("r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        if reader.fieldnames is None:
            raise ValueError(f"Missing CSV header: {csv_path}")
        missing = [col for col in REQUIRED_COLUMNS if col not in reader.fieldnames]
        if missing:
            raise ValueError(f"CSV schema mismatch in {csv_path}: missing {', '.join(missing)}")
        return rows


def f(value: str) -> float:
    try:
        return float(value)
    except Exception:
        return math.nan


def physical_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    out = []
    for row in rows:
        if int(float(row["is_solid"])) == 0:
            out.append(row)
    return out


def unique_sorted(values: list[float]) -> list[float]:
    return sorted(set(values))


def by_y(rows: list[dict[str, str]]) -> dict[float, list[dict[str, str]]]:
    grouped: dict[float, list[dict[str, str]]] = {}
    for row in rows:
        grouped.setdefault(f(row["y"]), []).append(row)
    return grouped


def pick_row_by_y(rows: list[dict[str, str]], target_y: float | None = None) -> list[dict[str, str]]:
    grouped = by_y(rows)
    if not grouped:
        return []
    ys = unique_sorted(list(grouped))
    if target_y is None:
        target_y = ys[len(ys) // 2]
    best_y = min(ys, key=lambda yy: abs(yy - target_y))
    return sorted(grouped[best_y], key=lambda r: f(r["x"]))


def mean(values: list[float]) -> float:
    values = [v for v in values if math.isfinite(v)]
    if not values:
        return math.nan
    return statistics.fmean(values)


def row_dx(profile: list[dict[str, str]]) -> float:
    xs = [f(row["x"]) for row in profile]
    diffs = [b - a for a, b in zip(xs, xs[1:]) if math.isfinite(a) and math.isfinite(b) and b > a]
    return min(diffs) if diffs else math.nan


def window_rows(
    profile: list[dict[str, str]],
    x0: float,
    x1: float,
    *,
    field: str,
) -> list[float]:
    return [f(row[field]) for row in profile if x0 <= f(row["x"]) <= x1 and math.isfinite(f(row[field]))]


def window_mean(profile: list[dict[str, str]], x0: float, x1: float, field: str) -> float:
    return mean(window_rows(profile, x0, x1, field=field))


def count_transition_cells(
    profile: list[dict[str, str]],
    *,
    x0: float,
    x1: float,
    field: str,
    low: float,
    high: float,
) -> int:
    count = 0
    for row in profile:
        x = f(row["x"])
        value = f(row[field])
        if x0 <= x <= x1 and math.isfinite(value) and low < value < high:
            count += 1
    return count


def gradient_peak_x(profile: list[dict[str, str]], field: str, x_min: float, x_max: float) -> float:
    best_x = math.nan
    best_score = -math.inf
    for idx in range(1, len(profile)):
        x0 = f(profile[idx - 1]["x"])
        x1 = f(profile[idx]["x"])
        if not (x_min <= x0 <= x_max and x_min <= x1 <= x_max):
            continue
        v0 = f(profile[idx - 1][field])
        v1 = f(profile[idx][field])
        dx = x1 - x0
        if dx <= 0.0:
            continue
        score = abs((v1 - v0) / dx)
        if score > best_score:
            best_score = score
            best_x = 0.5 * (x0 + x1)
    return best_x


def transition_width(profile: list[dict[str, str]], field: str) -> int:
    values = [f(row[field]) for row in profile]
    finite = [v for v in values if math.isfinite(v)]
    if not finite:
        return 0
    lo = min(finite)
    hi = max(finite)
    if not math.isfinite(lo) or not math.isfinite(hi) or hi <= lo:
        return 0
    low_band = lo + 0.05 * (hi - lo)
    high_band = hi - 0.05 * (hi - lo)
    return sum(1 for value in values if math.isfinite(value) and low_band <= value <= high_band)


def transition_is_monotone(profile: list[dict[str, str]], field: str) -> bool:
    values = [f(row[field]) for row in profile if math.isfinite(f(row[field]))]
    if len(values) < 2:
        return False
    return all(values[idx] <= values[idx + 1] + 1e-8 for idx in range(len(values) - 1))


def front_position_from_w(profile: list[dict[str, str]], threshold: float = 0.5) -> float:
    xs = [f(row["x"]) for row in profile]
    ws = [f(row["w"]) for row in profile]
    front = math.nan
    for x, w in zip(xs, ws):
        if math.isfinite(x) and math.isfinite(w) and w <= threshold:
            front = x
        else:
            break
    return front


def front_position_from_p(profile: list[dict[str, str]], threshold: float) -> float:
    xs = [f(row["x"]) for row in profile]
    ps = [f(row["p"]) for row in profile]
    for idx in range(len(ps) - 1, -1, -1):
        if math.isfinite(ps[idx]) and ps[idx] >= threshold:
            return xs[idx]
    return math.nan


def fit_line(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    if len(xs) != len(ys) or len(xs) < 2:
        return math.nan, math.nan, math.nan
    x_mean = statistics.fmean(xs)
    y_mean = statistics.fmean(ys)
    sxx = sum((x - x_mean) ** 2 for x in xs)
    sxy = sum((x - x_mean) * (y - y_mean) for x, y in zip(xs, ys))
    if sxx <= 0:
        return math.nan, math.nan, math.nan
    slope = sxy / sxx
    intercept = y_mean - slope * x_mean
    ss_tot = sum((y - y_mean) ** 2 for y in ys)
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(xs, ys))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0
    return slope, intercept, r2


def latest_csv(case_dir: Path) -> Path:
    final = case_dir / "final_results.csv"
    if final.exists():
        return final
    candidates = sorted(case_dir.glob("step_*_time_*.csv"), key=lambda path: frame_time(path))
    if not candidates:
        raise FileNotFoundError(f"No CSV outputs found in {case_dir}")
    return candidates[-1]


def frame_time(path: Path) -> float:
    stem = path.stem
    if stem == "final_results":
        return math.inf
    if stem == "step_0_initial":
        return 0.0
    if "_time_" in stem:
        try:
            return float(stem.rsplit("_time_", 1)[1])
        except Exception:
            return math.nan
    return math.nan


def load_case_dir(output_root: Path, case: str) -> Path:
    case_dir = output_root / case
    if not case_dir.exists():
        raise FileNotFoundError(f"Missing output directory: {case_dir}")
    return case_dir


def close_to(value: float, expected: float, rel_tol: float) -> bool:
    return math.isfinite(value) and abs(value - expected) <= rel_tol * max(1.0, abs(expected))


def sod_check(case_dir: Path) -> CaseResult:
    rows = physical_rows(load_rows(latest_csv(case_dir)))
    if not rows:
        return CaseResult("sod", False, "No physical rows found")

    profile = pick_row_by_y(rows)
    if not profile:
        return CaseResult("sod", False, "No 1D profile found")

    dx = row_dx(profile)
    if not math.isfinite(dx):
        return CaseResult("sod", False, "Could not infer grid spacing")

    contact_x = gradient_peak_x(profile, "rho", 0.60, 0.76)
    shock_x = gradient_peak_x(profile, "p", 0.78, 0.94)
    if not math.isfinite(contact_x) or not math.isfinite(shock_x):
        return CaseResult("sod", False, "Could not locate contact/shock")

    p_left = window_mean(profile, max(0.0, contact_x - 0.12), contact_x - 0.03, "p")
    rho_left = window_mean(profile, max(0.0, contact_x - 0.12), contact_x - 0.03, "rho")
    u_left = window_mean(profile, max(0.0, contact_x - 0.12), contact_x - 0.03, "u")
    p_right = window_mean(profile, contact_x + 0.03, shock_x - 0.03, "p")
    rho_right = window_mean(profile, contact_x + 0.03, shock_x - 0.03, "rho")
    u_right = window_mean(profile, contact_x + 0.03, shock_x - 0.03, "u")

    shock_width = count_transition_cells(
        profile,
        x0=shock_x - 0.05,
        x1=shock_x + 0.05,
        field="rho",
        low=SOD_EXPECTED["rho_far_right"],
        high=SOD_EXPECTED["rho_right"],
    )
    contact_width = count_transition_cells(
        profile,
        x0=contact_x - 0.06,
        x1=contact_x + 0.06,
        field="rho",
        low=SOD_EXPECTED["rho_right"],
        high=SOD_EXPECTED["rho_left"],
    )
    max_w = max(abs(f(row["w"])) for row in profile)

    checks = [
        (max_w <= 1e-9, f"W should remain zero, got max|W|={max_w:.3e}"),
        (close_to(p_left, SOD_EXPECTED["p_left"], 0.05), f"left plateau pressure out of range: {p_left:.4f}"),
        (close_to(rho_left, SOD_EXPECTED["rho_left"], 0.05), f"left plateau density out of range: {rho_left:.4f}"),
        (close_to(u_left, SOD_EXPECTED["u_left"], 0.05), f"left plateau velocity out of range: {u_left:.4f}"),
        (close_to(p_right, SOD_EXPECTED["p_right"], 0.05), f"right plateau pressure out of range: {p_right:.4f}"),
        (close_to(rho_right, SOD_EXPECTED["rho_right"], 0.05), f"right plateau density out of range: {rho_right:.4f}"),
        (close_to(u_right, SOD_EXPECTED["u_right"], 0.05), f"right plateau velocity out of range: {u_right:.4f}"),
        (abs(contact_x - SOD_EXPECTED["contact_x"]) <= 3.0 * dx, f"contact position out of range: {contact_x:.4f}"),
        (abs(shock_x - SOD_EXPECTED["shock_x"]) <= 3.0 * dx, f"shock position out of range: {shock_x:.4f}"),
        (5 <= contact_width <= 8, f"contact smearing should be 5-8 cells, got {contact_width}"),
        (3 <= shock_width <= 5, f"shock smearing should be 3-5 cells, got {shock_width}"),
    ]
    failed = [msg for ok, msg in checks if not ok]
    if failed:
        return CaseResult("sod", False, "; ".join(failed))

    return CaseResult(
        "sod",
        True,
        "p_left={:.4f}, rho_left={:.4f}, u_left={:.4f}, p_right={:.4f}, rho_right={:.4f}, u_right={:.4f}, "
        "contact_x={:.4f}, shock_x={:.4f}, contact_width={}, shock_width={}".format(
            p_left,
            rho_left,
            u_left,
            p_right,
            rho_right,
            u_right,
            contact_x,
            shock_x,
            contact_width,
            shock_width,
        ),
    )


def sorted_frames(case_dir: Path) -> list[Path]:
    csvs = list(case_dir.glob("step_*_time_*.csv"))
    initial = case_dir / "step_0_initial.csv"
    if initial.exists():
        csvs.append(initial)
    csvs.sort(key=frame_time)
    final = case_dir / "final_results.csv"
    if final.exists():
        csvs.append(final)
    return csvs


def cj_check(case_dir: Path) -> CaseResult:
    csvs = sorted_frames(case_dir)
    if not csvs:
        return CaseResult("cj", False, "No frames found")

    times: list[float] = []
    fronts: list[float] = []
    initial_profile: list[dict[str, str]] | None = None
    final_profile: list[dict[str, str]] | None = None

    for csv_path in csvs:
        rows = physical_rows(load_rows(csv_path))
        profile = pick_row_by_y(rows)
        if not profile:
            continue
        if initial_profile is None:
            initial_profile = profile
        final_profile = profile
        front = front_position_from_w(profile)
        if math.isfinite(front):
            times.append(frame_time(csv_path))
            fronts.append(front)

    if len(fronts) < 3:
        return CaseResult("cj", False, "Not enough usable frames to estimate the front speed")
    if initial_profile is None or final_profile is None:
        return CaseResult("cj", False, "Missing final profile")

    finite_pairs = [(t, x) for t, x in zip(times, fronts) if math.isfinite(t) and math.isfinite(x)]
    if len(finite_pairs) < 3:
        return CaseResult("cj", False, "Not enough finite front positions")

    times = [t for t, _ in finite_pairs]
    fronts = [x for _, x in finite_pairs]
    slope, intercept, r2 = fit_line(times, fronts)
    if not math.isfinite(slope):
        return CaseResult("cj", False, "Front-speed regression failed")

    d_cj = CJ_EXPECTED["d_cj"]
    speed_ok = abs(slope - d_cj) <= 0.02 * d_cj
    linear_ok = r2 >= 0.98
    motion_ok = all(fronts[idx] <= fronts[idx + 1] + 1e-9 for idx in range(len(fronts) - 1))

    dx = row_dx(final_profile)
    zone_width = transition_width(final_profile, "w")
    final_front = front_position_from_w(final_profile)
    behind_w = window_mean(final_profile, max(0.0, final_front - 0.30), max(0.0, final_front - 0.10), "w")
    ahead_w = window_mean(final_profile, final_front + 0.10, min(f(final_profile[-1]["x"]), final_front + 0.30), "w")
    behind_p = window_mean(final_profile, max(0.0, final_front - 0.30), max(0.0, final_front - 0.10), "p")
    behind_pressure_values = [
        f(row["p"]) for row in final_profile if math.isfinite(f(row["p"])) and f(row["x"]) <= final_front - 0.10
    ]
    if behind_pressure_values:
        p_std = statistics.pstdev(behind_pressure_values)
        p_mean = mean(behind_pressure_values)
        oscillation_ok = math.isfinite(p_mean) and (p_std / max(abs(p_mean), 1e-12) <= 0.10)
    else:
        p_std = math.nan
        p_mean = math.nan
        oscillation_ok = False

    detonator_p = window_mean(initial_profile, 0.0, 0.03, "p")
    detonator_rho = window_mean(initial_profile, 0.0, 0.03, "rho")
    detonator_w = window_mean(initial_profile, 0.0, 0.03, "w")
    unreacted_p = window_mean(initial_profile, 0.04, 0.20, "p")
    unreacted_rho = window_mean(initial_profile, 0.04, 0.20, "rho")
    unreacted_w = window_mean(initial_profile, 0.04, 0.20, "w")

    checks = [
        (close_to(detonator_p, CJ_EXPECTED["p_cj"], 0.05), f"initial detonator pressure out of range: {detonator_p:.4f}"),
        (close_to(detonator_rho, CJ_EXPECTED["detonator_rho"], 0.05), f"initial detonator density out of range: {detonator_rho:.4f}"),
        (math.isfinite(detonator_w) and detonator_w <= 0.05, f"initial detonator W should be 0, got {detonator_w:.4f}"),
        (abs(unreacted_p) <= 0.05 * max(1.0, CJ_EXPECTED["p_cj"]), f"initial unreacted pressure should be near zero, got {unreacted_p:.4f}"),
        (close_to(unreacted_rho, CJ_EXPECTED["rho0"], 0.05), f"initial unreacted density out of range: {unreacted_rho:.4f}"),
        (math.isfinite(unreacted_w) and unreacted_w >= 0.95, f"initial unreacted W should be 1, got {unreacted_w:.4f}"),
        (speed_ok, f"front speed out of range: {slope:.4f}"),
        (linear_ok, f"front motion is too nonlinear: r2={r2:.4f}"),
        (motion_ok, f"front should not retreat: {fronts}"),
        (zone_width >= 5 and zone_width <= 10, f"reaction zone width should be 5-10 cells, got {zone_width}"),
        (math.isfinite(behind_w) and behind_w <= 0.05, f"W behind the front should go to zero, got {behind_w:.4f}"),
        (math.isfinite(ahead_w) and ahead_w >= 0.90, f"W ahead of the front should remain near one, got {ahead_w:.4f}"),
        (oscillation_ok, f"pressure oscillations behind the front are too large: mean={p_mean:.4f}, std={p_std:.4f}"),
    ]
    failed = [msg for ok, msg in checks if not ok]
    if failed:
        return CaseResult("cj", False, "; ".join(failed))

    return CaseResult(
        "cj",
        True,
        "speed={:.4f}, r2={:.4f}, front_final={:.4f}, zone_width={}, W_behind={:.4f}, W_ahead={:.4f}, "
        "p_behind_mean={:.4f}".format(
            slope,
            r2,
            final_front,
            zone_width,
            behind_w,
            ahead_w,
            behind_p,
        ),
    )


def finite_numeric(row: dict[str, str]) -> bool:
    for key in ("rho", "u", "v", "p", "w", "temperature"):
        value = f(row[key])
        if not math.isfinite(value):
            return False
    return True


def corner_check(case_dir: Path) -> CaseResult:
    csvs = sorted_frames(case_dir)
    if not csvs:
        return CaseResult("corner", False, "No frames found")

    last_rows = physical_rows(load_rows(csvs[-1]))
    if not last_rows:
        return CaseResult("corner", False, "No physical rows in final frame")

    if any(not finite_numeric(row) for row in last_rows):
        return CaseResult("corner", False, "Non-finite physical values found")

    rho_min = min(f(r["rho"]) for r in last_rows)
    p_min = min(f(r["p"]) for r in last_rows)
    if rho_min <= 0.0 or p_min <= 0.0:
        return CaseResult("corner", False, f"Negative or zero state detected: rho_min={rho_min:.4f}, p_min={p_min:.4f}")

    by_level = by_y(last_rows)
    y_levels = sorted(by_level)
    if len(y_levels) < 10:
        return CaseResult("corner", False, "Too few y-levels for curvature analysis")

    channel_profile = pick_row_by_y(last_rows, CORNER_EXPECTED["channel_top_y"] + 0.5)
    if not channel_profile:
        return CaseResult("corner", False, "Could not identify the channel profile")

    first_rows = physical_rows(load_rows(csvs[0]))
    first_profile = pick_row_by_y(first_rows, CORNER_EXPECTED["channel_top_y"] + 0.5)
    if not first_profile:
        return CaseResult("corner", False, "Could not identify the initial channel profile")

    channel_frames = []
    for csv_path in csvs:
        rows = physical_rows(load_rows(csv_path))
        profile = pick_row_by_y(rows, CORNER_EXPECTED["channel_top_y"] + 0.5)
        if not profile:
            continue
        front = front_position_from_w(profile)
        if math.isfinite(front) and front <= CORNER_EXPECTED["channel_exit_x"] - 0.05:
            channel_frames.append((frame_time(csv_path), front))

    if len(channel_frames) < 3:
        return CaseResult("corner", False, "Not enough in-channel frames to estimate speed")

    channel_times = [t for t, _ in channel_frames]
    channel_fronts = [x for _, x in channel_frames]
    finite_pairs = [(t, x) for t, x in zip(channel_times, channel_fronts) if math.isfinite(t) and math.isfinite(x)]
    if len(finite_pairs) < 3:
        return CaseResult("corner", False, "Not enough finite in-channel frames to estimate speed")
    channel_times = [t for t, _ in finite_pairs]
    channel_fronts = [x for _, x in finite_pairs]
    slope, intercept, r2 = fit_line(channel_times, channel_fronts)
    if not math.isfinite(slope):
        return CaseResult("corner", False, "Channel front regression failed")

    d_cj = CORNER_EXPECTED["d_cj"]
    speed_ok = abs(slope - d_cj) <= 0.02 * d_cj
    linear_ok = r2 >= 0.98

    final_by_y = by_level
    front_positions: list[float] = []
    for y in y_levels:
        profile = sorted(final_by_y[y], key=lambda r: f(r["x"]))
        pos = front_position_from_w(profile)
        if math.isfinite(pos):
            front_positions.append(pos)

    if len(front_positions) < 10:
        return CaseResult("corner", False, "Could not identify enough front positions")

    curvature = max(front_positions) - min(front_positions)
    curved_ok = curvature >= 0.20

    top_window = [
        row
        for row in last_rows
        if f(row["x"]) >= 4.5 and f(row["y"]) >= 3.15
    ]
    bottom_window = [
        row
        for row in last_rows
        if f(row["x"]) >= 4.5 and f(row["y"]) <= 1.0
    ]
    p_top = mean([f(row["p"]) for row in top_window])
    p_bottom = mean([f(row["p"]) for row in bottom_window])
    w_top = mean([f(row["w"]) for row in top_window])
    w_bottom = mean([f(row["w"]) for row in bottom_window])
    initial_p = window_mean(first_profile, 0.0, 0.03, "p")
    initial_rho = window_mean(first_profile, 0.0, 0.03, "rho")
    initial_w = window_mean(first_profile, 0.0, 0.03, "w")
    dead_zone_ok = (
        math.isfinite(p_top)
        and math.isfinite(p_bottom)
        and math.isfinite(w_top)
        and math.isfinite(w_bottom)
        and p_top > p_bottom * 1.2
        and w_bottom > w_top + 0.1
    )

    checks = [
        (close_to(initial_p, CORNER_EXPECTED["p_cj"], 0.05), f"initial detonator pressure out of range: {initial_p:.4f}"),
        (close_to(initial_rho, CJ_EXPECTED["detonator_rho"], 0.05), f"initial detonator density out of range: {initial_rho:.4f}"),
        (math.isfinite(initial_w) and initial_w <= 0.05, f"initial detonator W should be 0, got {initial_w:.4f}"),
        (speed_ok, f"in-channel front speed out of range: {slope:.4f}"),
        (linear_ok, f"in-channel front motion is too nonlinear: r2={r2:.4f}"),
        (curved_ok, f"front curvature too small: {curvature:.4f}"),
        (dead_zone_ok, f"dead-zone contrast too weak: p_top={p_top:.4f}, p_bottom={p_bottom:.4f}, w_top={w_top:.4f}, w_bottom={w_bottom:.4f}"),
    ]
    failed = [msg for ok, msg in checks if not ok]
    if failed:
        return CaseResult("corner", False, "; ".join(failed))

    return CaseResult(
        "corner",
        True,
        "speed={:.4f}, r2={:.4f}, curvature={:.4f}, p_top={:.4f}, p_bottom={:.4f}, w_top={:.4f}, w_bottom={:.4f}".format(
            slope,
            r2,
            curvature,
            p_top,
            p_bottom,
            w_top,
            w_bottom,
        ),
    )


def verify_case(case: str, output_root: Path) -> CaseResult:
    case_dir = load_case_dir(output_root, case)
    if case == "sod":
        return sod_check(case_dir)
    if case == "cj":
        return cj_check(case_dir)
    if case == "corner":
        return corner_check(case_dir)
    raise ValueError(f"Unknown case: {case}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify Mader 2DE outputs")
    parser.add_argument("--case", choices=["sod", "cj", "corner", "all"], required=True)
    parser.add_argument("--output-root", type=str, default="tests/mader_2de/output")
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[2]
    output_root = (root / args.output_root).resolve()
    cases = ["sod", "cj", "corner"] if args.case == "all" else [args.case]

    results = [verify_case(case, output_root) for case in cases]
    for result in results:
        status = "PASS" if result.passed else "FAIL"
        print(f"{status} {result.name}: {result.message}")
    if any(not result.passed for result in results):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
