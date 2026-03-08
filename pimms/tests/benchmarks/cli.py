from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import shlex
import shutil
import statistics
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass
class BenchmarkRun:
    scenario: int
    run_index: int
    started_at: str
    finished_at: str
    seconds: float
    returncode: int
    success: bool
    final_energy: str | None
    timed_out: bool


@dataclass
class ScenarioSummary:
    scenario: int
    attempted_runs: int
    successful_runs: int
    failed_runs: int
    min_seconds: float | None
    mean_seconds: float | None
    max_seconds: float | None
    stdev_seconds: float | None


def _repo_root() -> Path:
    current = Path(__file__).resolve()
    for parent in current.parents:
        if (parent / "setup.py").exists() and (parent / "scripts" / "PIMMS").exists():
            return parent
    raise RuntimeError(f"Unable to locate repository root from {current}")


def _simulation_tests_root() -> Path:
    return Path(__file__).resolve().parents[1] / "simulation_tests"


def _benchmarks_root() -> Path:
    return Path(__file__).resolve().parents[1] / "benchmarks"


def _resolve_pimms_command(user_command: str | None) -> list[str]:
    if user_command:
        return shlex.split(user_command)

    pimms_exe = shutil.which("PIMMS")
    if pimms_exe:
        return [pimms_exe]

    script_path = _repo_root() / "scripts" / "PIMMS"
    return [sys.executable, str(script_path)]


def _discover_scenarios(sim_root: Path) -> list[int]:
    scenarios: list[int] = []
    for path in sim_root.glob("test_*"):
        if not path.is_dir():
            continue
        suffix = path.name.removeprefix("test_")
        if suffix.isdigit():
            scenarios.append(int(suffix))
    return sorted(scenarios)


def _parse_scenario_selection(value: str, available: Iterable[int]) -> list[int]:
    available_set = set(available)
    selected: set[int] = set()

    for token in (chunk.strip() for chunk in value.split(",") if chunk.strip()):
        if "-" in token:
            start_s, end_s = token.split("-", 1)
            start_i = int(start_s)
            end_i = int(end_s)
            if end_i < start_i:
                raise ValueError(f"Invalid scenario range: {token}")
            for number in range(start_i, end_i + 1):
                selected.add(number)
            continue

        selected.add(int(token))

    missing = sorted(s for s in selected if s not in available_set)
    if missing:
        raise ValueError(f"Unknown scenario(s): {missing}; available scenarios: {sorted(available_set)}")

    return sorted(selected)


def _cleanup_generated_outputs(test_dir: Path) -> None:
    for path in test_dir.iterdir():
        if not path.is_file():
            continue

        if path.suffix in {".prm", ".kf"}:
            continue

        if path.name.startswith("KEYFILE"):
            continue

        if path.suffix in {".dat", ".xtc", ".pdb", ".txt", ".lat"} or path.name == "restart.pimms":
            path.unlink(missing_ok=True)


def _extract_final_energy(test_dir: Path) -> str | None:
    energy_path = test_dir / "ENERGY.dat"
    if not energy_path.exists():
        return None

    lines = [ln.strip() for ln in energy_path.read_text().splitlines() if ln.strip()]
    if not lines:
        return None

    fields = lines[-1].split()
    if len(fields) < 2:
        return None

    return fields[1]


def _run_single(
    cmd_base: list[str],
    scenario: int,
    run_index: int,
    timeout_seconds: int,
    keep_outputs: bool,
    sim_root: Path,
    log_dir: Path,
) -> BenchmarkRun:
    scenario_dir = sim_root / f"test_{scenario}"

    if not keep_outputs:
        _cleanup_generated_outputs(scenario_dir)

    command = cmd_base + ["-k", "KEYFILE.kf"]
    started_at = dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")
    started = time.perf_counter()
    timed_out = False
    try:
        result = subprocess.run(
            command,
            cwd=str(scenario_dir),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=timeout_seconds,
            check=False,
        )
        returncode = result.returncode
    except subprocess.TimeoutExpired as exc:
        # Mark timeout as a failed run and continue the benchmark suite.
        timed_out = True
        _ = exc
        returncode = 124

    elapsed = time.perf_counter() - started
    finished_at = dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")

    final_energy = _extract_final_energy(scenario_dir)

    return BenchmarkRun(
        scenario=scenario,
        run_index=run_index,
        started_at=started_at,
        finished_at=finished_at,
        seconds=elapsed,
        returncode=returncode,
        success=returncode == 0,
        final_energy=final_energy,
        timed_out=timed_out,
    )


def _initialize_status_log(path: Path, cmd_base: list[str], scenarios: list[int], repeats: int, timeout: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    started_at = dt.datetime.now(dt.timezone.utc).astimezone().isoformat(timespec="seconds")
    with path.open("w") as handle:
        handle.write("PIMMS benchmark status log\n")
        handle.write(f"started_at: {started_at}\n")
        handle.write(f"command: {' '.join(cmd_base)}\n")
        handle.write(f"scenarios: {scenarios}\n")
        handle.write(f"repeats: {repeats}\n")
        handle.write(f"timeout_seconds: {timeout}\n")
        handle.write("-" * 72 + "\n")


def _append_status_log(path: Path, run: BenchmarkRun) -> None:
    status = "ok" if run.success else "fail"
    timeout_note = " timeout" if run.timed_out else ""
    energy = run.final_energy if run.final_energy is not None else "n/a"
    with path.open("a") as handle:
        handle.write(
            f"{run.finished_at} | scenario={run.scenario} run={run.run_index} |"
            f" status={status}{timeout_note} returncode={run.returncode}"
            f" elapsed_s={run.seconds:.3f} final_energy={energy}"
            f" started_at={run.started_at}\n"
        )


def _summarize(runs: list[BenchmarkRun]) -> list[ScenarioSummary]:
    by_scenario: dict[int, list[BenchmarkRun]] = {}
    for run in runs:
        by_scenario.setdefault(run.scenario, []).append(run)

    summaries: list[ScenarioSummary] = []
    for scenario in sorted(by_scenario):
        scenario_runs = by_scenario[scenario]
        success_times = [r.seconds for r in scenario_runs if r.success]
        failed = sum(1 for r in scenario_runs if not r.success)
        summaries.append(
            ScenarioSummary(
                scenario=scenario,
                attempted_runs=len(scenario_runs),
                successful_runs=len(success_times),
                failed_runs=failed,
                min_seconds=min(success_times) if success_times else None,
                mean_seconds=statistics.mean(success_times) if success_times else None,
                max_seconds=max(success_times) if success_times else None,
                stdev_seconds=statistics.stdev(success_times) if len(success_times) > 1 else None,
            )
        )

    return summaries


def _print_run_table(runs: list[BenchmarkRun]) -> None:
    print("\nPer-run results")
    print("scenario run  status  seconds   returncode final_energy")
    for run in runs:
        status = "ok" if run.success else "fail"
        if run.timed_out:
            status = "timeout"
        energy = run.final_energy if run.final_energy is not None else "n/a"
        print(
            f"{run.scenario:>8} {run.run_index:>3}  {status:<6} {run.seconds:>8.3f}"
            f" {run.returncode:>10} {energy}"
        )


def _print_summary_table(summaries: list[ScenarioSummary]) -> None:
    print("\nScenario summary")
    print("scenario runs ok fail min_s mean_s max_s stdev_s")
    for item in summaries:
        min_s = f"{item.min_seconds:.3f}" if item.min_seconds is not None else "n/a"
        mean_s = f"{item.mean_seconds:.3f}" if item.mean_seconds is not None else "n/a"
        max_s = f"{item.max_seconds:.3f}" if item.max_seconds is not None else "n/a"
        stdev_s = f"{item.stdev_seconds:.3f}" if item.stdev_seconds is not None else "n/a"
        print(
            f"{item.scenario:>8} {item.attempted_runs:>4} {item.successful_runs:>2}"
            f" {item.failed_runs:>4} {min_s:>5} {mean_s:>6} {max_s:>5} {stdev_s:>7}"
        )


def _write_csv(path: Path, runs: list[BenchmarkRun]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "scenario",
                "run",
                "started_at",
                "finished_at",
                "status",
                "seconds",
                "returncode",
                "final_energy",
                "timed_out",
            ],
        )
        writer.writeheader()
        for run in runs:
            writer.writerow(
                {
                    "scenario": run.scenario,
                    "run": run.run_index,
                    "started_at": run.started_at,
                    "finished_at": run.finished_at,
                    "status": "ok" if run.success else "fail",
                    "seconds": f"{run.seconds:.6f}",
                    "returncode": run.returncode,
                    "final_energy": run.final_energy or "",
                    "timed_out": str(run.timed_out).lower(),
                }
            )


def _write_json(path: Path, runs: list[BenchmarkRun], summaries: list[ScenarioSummary]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    payload = {
        "runs": [
            {
                "scenario": run.scenario,
                "run": run.run_index,
                "started_at": run.started_at,
                "finished_at": run.finished_at,
                "status": "ok" if run.success else "fail",
                "seconds": run.seconds,
                "returncode": run.returncode,
                "final_energy": run.final_energy,
                "timed_out": run.timed_out,
            }
            for run in runs
        ],
        "summary": [
            {
                "scenario": item.scenario,
                "attempted_runs": item.attempted_runs,
                "successful_runs": item.successful_runs,
                "failed_runs": item.failed_runs,
                "min_seconds": item.min_seconds,
                "mean_seconds": item.mean_seconds,
                "max_seconds": item.max_seconds,
                "stdev_seconds": item.stdev_seconds,
            }
            for item in summaries
        ],
    }

    path.write_text(json.dumps(payload, indent=2))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Benchmark full PIMMS CLI simulations across selected scenario fixtures in "
            "pimms/tests/simulation_tests."
        )
    )
    parser.add_argument(
        "--scenarios",
        default="all",
        help=(
            "Scenario selector: 'all', comma lists (e.g. 1,2,4), or ranges "
            "(e.g. 1-4,8,10-12)."
        ),
    )
    parser.add_argument("--repeats", type=int, default=1, help="Runs per selected scenario.")
    parser.add_argument("--timeout", type=int, default=1800, help="Per-run timeout in seconds.")
    parser.add_argument(
        "--pimms-command",
        default=None,
        help=(
            "Explicit command to run PIMMS, for example: 'PIMMS' or "
            "'/usr/bin/python scripts/PIMMS'."
        ),
    )
    parser.add_argument(
        "--keep-outputs",
        action="store_true",
        help="Do not remove generated simulation outputs between runs.",
    )
    parser.add_argument(
        "--log-dir",
        default=None,
        help="Directory for benchmark logs (default: pimms/tests/benchmarks/benchmark_logs).",
    )
    parser.add_argument(
        "--status-log",
        default=None,
        help=(
            "Path to a single run-status logfile. Default: <log-dir>/benchmark_status_<YYYYmmdd_HHMMSS>.log. "
            "This file includes per-run timestamp + pass/fail status."
        ),
    )
    parser.add_argument("--csv-out", default=None, help="Optional CSV output path for per-run rows.")
    parser.add_argument("--json-out", default=None, help="Optional JSON output path for full report.")
    parser.add_argument(
        "--stop-on-error",
        action="store_true",
        help="Stop immediately if any scenario run fails.",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.repeats < 1:
        parser.error("--repeats must be >= 1")
    if args.timeout < 1:
        parser.error("--timeout must be >= 1")

    sim_root = _simulation_tests_root()
    available = _discover_scenarios(sim_root)
    if not available:
        print(f"No scenario directories found in {sim_root}", file=sys.stderr)
        return 2

    if args.scenarios.strip().lower() == "all":
        selected_scenarios = available
    else:
        try:
            selected_scenarios = _parse_scenario_selection(args.scenarios, available)
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            return 2

    cmd_base = _resolve_pimms_command(args.pimms_command)
    run_stamp = dt.datetime.now(dt.timezone.utc).astimezone().strftime("%Y%m%d_%H%M%S")
    default_log_dir = _benchmarks_root() / "benchmark_logs"
    log_dir = Path(args.log_dir) if args.log_dir else default_log_dir
    log_dir.mkdir(parents=True, exist_ok=True)
    status_log_path = (
        Path(args.status_log)
        if args.status_log
        else (log_dir / f"benchmark_status_{run_stamp}.log")
    )
    _initialize_status_log(
        path=status_log_path,
        cmd_base=cmd_base,
        scenarios=selected_scenarios,
        repeats=args.repeats,
        timeout=args.timeout,
    )

    print(f"PIMMS command: {' '.join(cmd_base)}")
    print(f"Scenarios: {selected_scenarios}")
    print(f"Repeats per scenario: {args.repeats}")
    print(f"Per-run timeout: {args.timeout} s")
    print(f"Log directory: {log_dir}")
    print(f"Status log: {status_log_path}")

    runs: list[BenchmarkRun] = []
    for scenario in selected_scenarios:
        for run_index in range(1, args.repeats + 1):
            print(f"Running scenario {scenario} (run {run_index}/{args.repeats})...")
            run = _run_single(
                cmd_base=cmd_base,
                scenario=scenario,
                run_index=run_index,
                timeout_seconds=args.timeout,
                keep_outputs=args.keep_outputs,
                sim_root=sim_root,
                log_dir=log_dir,
            )
            runs.append(run)

            state = "ok" if run.success else "fail"
            if run.timed_out:
                state = "timeout"
            energy = run.final_energy if run.final_energy is not None else "n/a"
            print(
                f"  -> {state}: {run.seconds:.3f}s, returncode={run.returncode}, "
                f"final_energy={energy}"
            )
            _append_status_log(status_log_path, run)

            if args.stop_on_error and not run.success:
                print("Stopping early due to --stop-on-error")
                summaries = _summarize(runs)
                _print_run_table(runs)
                _print_summary_table(summaries)
                return 1

    summaries = _summarize(runs)
    _print_run_table(runs)
    _print_summary_table(summaries)

    if args.csv_out:
        csv_path = Path(args.csv_out)
        _write_csv(csv_path, runs)
        print(f"Wrote CSV report: {csv_path}")

    if args.json_out:
        json_path = Path(args.json_out)
        _write_json(json_path, runs, summaries)
        print(f"Wrote JSON report: {json_path}")

    failures = sum(1 for run in runs if not run.success)
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
