# benchmark

Command-line benchmarking helper for running full PIMMS CLI simulations over the
existing `pimms/tests/simulation_tests/test_*` scenarios.

Simulation stdout/stderr is discarded during benchmark execution. Instead, a
single status logfile is written with per-run timestamps and pass/fail status.

## Run

```bash
# from the benchmarking directory (recommended)
python cli.py

# from root
python -m pimms.tests.benchmarks --scenarios all --repeats 3
```

## Useful options

- `--scenarios`: `all`, comma lists (`1,2,4`) and ranges (`1-4,8,10-12`)
- `--repeats`: number of runs per scenario
- `--timeout`: per-run timeout in seconds
- `--pimms-command`: explicit executable string if `PIMMS` is not on `PATH`
- `--log-dir`: directory for benchmark logs
- `--status-log`: explicit path for the single status logfile
- `--csv-out`: write per-run results as CSV
- `--json-out`: write full report with per-run data + scenario summaries
- `--stop-on-error`: stop immediately when one run fails
- `--keep-outputs`: keep generated files between benchmark runs

Default status logfile:

- default directory: `pimms/tests/benchmarks/benchmark_logs/`
- default filename: `benchmark_status_YYYYmmdd_HHMMSS.log`

## Example

```bash
# from this directory
python cli.py

python -m pimms.tests.benchmarks \
  --scenarios 1-6,10 \
  --repeats 5 \
  --csv-out /tmp/pimms_benchmark.csv \
  --json-out /tmp/pimms_benchmark.json
```
