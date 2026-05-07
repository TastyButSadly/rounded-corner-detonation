# Rounded-Corner Study

Run from repository root:

```bash
python3 tests/mader_2de/corner_fast_study.py \
  --solver build/GodunovSolver \
  --base-config tests/mader_2de/configs/corner_radius_study_base.ini \
  --output-root tests/mader_2de/output_corner_radius_run \
  --radii 0,0.25,0.5,1,1.5,2 \
  --nx 560 --ny 320 \
  --dt-out 0.2
```

Analyze existing outputs only:

```bash
python3 tests/mader_2de/corner_fast_study.py \
  --output-root tests/mader_2de/output_corner_radius_run \
  --radii 0,0.25,0.5,1,1.5,2 \
  --nx 560 --ny 320 \
  --skip-run
```

Edit config:

```text
tests/mader_2de/configs/corner_radius_study_base.ini
```

Main keys: `Nx`, `Ny`, `tmax`, `dt_out`, `step_x_end`, `step_y_end`,
`corner_radius`, `reaction_rate`, `reaction_activation_energy`,
`reaction_heat_release`.

Outputs: `summary.csv`, `critical_radius.txt`, `plots/S_dead_vs_R.png`,
`plots/S_hat_vs_R.png`, `R.../plots/*_dead_zone.png`.
