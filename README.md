# Rounded-Corner Detonation Diffraction

Numerical study of detonation diffraction through a channel turn with a rounded inner corner.

## Build

```bash
cmake -S . -B build
cmake --build build
```

The executable is `build/GodunovSolver`.

## Run Radius Study

```bash
python3 tests/mader_2de/corner_fast_study.py \
  --solver build/GodunovSolver \
  --base-config tests/mader_2de/configs/corner_radius_study_base.ini \
  --output-root tests/mader_2de/output_corner_radius_run \
  --radii 0,0.25,0.5,1,1.5,2 \
  --nx 560 --ny 320 \
  --dt-out 0.2
```

Edit the base config at:

```text
tests/mader_2de/configs/corner_radius_study_base.ini
```

Main outputs are written under `--output-root`:

- `summary.csv`
- `critical_radius.txt`
- `plots/S_dead_vs_R.png`
- `plots/S_hat_vs_R.png`
