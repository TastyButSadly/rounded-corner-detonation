#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class CaseSpec:
    name: str
    config: Path


CASES = [
    CaseSpec("sod", Path("configs/sod_2de.ini")),
    CaseSpec("cj", Path("configs/cj_1d.ini")),
    CaseSpec("corner", Path("configs/corner_turning.ini")),
]


def repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def default_solver_candidates(root: Path) -> list[Path]:
    return [
        root / "build" / "GodunovSolver",
        root / "build-codex" / "GodunovSolver",
        root / "build-mpi" / "GodunovSolver",
    ]


def run_cmd(cmd: list[str], *, cwd: Path, env: dict[str, str]) -> int:
    proc = subprocess.run(cmd, cwd=str(cwd), env=env)
    return proc.returncode


def resolve_solver(root: Path, override: str | None) -> Path:
    if override:
        solver = Path(override).expanduser()
        if not solver.is_absolute():
            solver = (root / solver).resolve()
        if not solver.exists():
            raise SystemExit(f"Solver not found: {solver}")
        return solver

    for candidate in default_solver_candidates(root):
        if candidate.exists():
            return candidate
    candidates = ", ".join(str(p) for p in default_solver_candidates(root))
    raise SystemExit(f"Solver not found. Tried: {candidates}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Run all Mader 2DE tests")
    parser.add_argument("--solver", type=str, default=None, help="Path to GodunovSolver binary")
    parser.add_argument("--output-root", type=str, default="tests/mader_2de/output")
    parser.add_argument("--cases", type=str, default="sod,cj,corner")
    parser.add_argument("--keep-existing", action="store_true")
    args = parser.parse_args()

    root = repo_root()
    solver = resolve_solver(root, args.solver)
    output_root = (root / args.output_root).resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    selected = {part.strip() for part in args.cases.split(",") if part.strip()}
    case_map = {case.name: case for case in CASES}
    unknown = sorted(selected - set(case_map))
    if unknown:
        raise SystemExit(f"Unknown cases requested: {', '.join(unknown)}")

    env = os.environ.copy()
    env.setdefault("MPLBACKEND", "Agg")
    env.setdefault("MPLCONFIGDIR", str((root / "tests" / "mader_2de" / ".mplcache").resolve()))

    failures: list[str] = []

    for case_name in [case.name for case in CASES if case.name in selected]:
        case = case_map[case_name]
        config = (root / "tests" / "mader_2de" / case.config).resolve()
        case_output = output_root / case.name
        if case_output.exists() and not args.keep_existing:
            shutil.rmtree(case_output)
        case_output.mkdir(parents=True, exist_ok=True)

        print(f"[run] {case.name}: {solver} {config} {case_output}")
        solver_code = run_cmd([str(solver), str(config), str(case_output)], cwd=root, env=env)
        if solver_code != 0:
            failures.append(f"{case.name}: solver exited with code {solver_code}")
            continue

        print(f"[plot] {case.name}")
        plot_code = run_cmd(
            [
                sys.executable,
                str((root / "tests" / "mader_2de" / "plot_results.py").resolve()),
                "--case",
                case.name,
                "--output-root",
                str(output_root),
            ],
            cwd=root,
            env=env,
        )
        if plot_code != 0:
            failures.append(f"{case.name}: plotting exited with code {plot_code}")

        print(f"[verify] {case.name}")
        verify_code = run_cmd(
            [
                sys.executable,
                str((root / "tests" / "mader_2de" / "verify.py").resolve()),
                "--case",
                case.name,
                "--output-root",
                str(output_root),
            ],
            cwd=root,
            env=env,
        )
        if verify_code != 0:
            failures.append(f"{case.name}: verification failed")

    if failures:
        print("[summary] completed with failures:")
        for failure in failures:
            print(f"  - {failure}")
        return 1

    print("[summary] all requested cases completed successfully")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
