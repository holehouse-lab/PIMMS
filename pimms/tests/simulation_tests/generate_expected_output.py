from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path


DEFAULT_FILES_TO_CAPTURE = ("ENERGY.dat", 
							"CLUSTERS.dat", 
							"CHAIN_0_CLUSTERS.dat", 
							"CHAIN_0_DISTANCE_MAP.dat",
                            "CHAIN_0_INTSCAL_SQUARED.dat",
                            "CHAIN_0_LR_CLUSTERS.dat",
                            "CHAIN_0_SCALING_INFORMATION.dat",
                            "CLUSTER_AREA.dat",
                            "CLUSTER_ASPH.dat",
                            "CLUSTER_DEN.dat",
                            "CLUSTER_RADIAL_DENSITY_PROFILE.dat",
                            "CLUSTER_RG.dat",
                            "CLUSTER_VOL.dat",
                            "END_TO_END_DIST.dat",
                            "LR_CLUSTER_AREA.dat",
                            "LR_CLUSTER_ASPH.dat",
                            "LR_CLUSTER_DEN.dat",
                            "LR_CLUSTER_RADIAL_DENSITY_PROFILE.dat",
                            "LR_CLUSTER_RG.dat",
                            "LR_CLUSTER_VOL.dat",
                            "LR_CLUSTERS.dat",
                            "MOVE_FREQS.dat",
                            "NUM_CLUSTERS.dat",
                            "NUM_LR_CLUSTERS.dat",
                            "RES_TO_RES_DIST.dat",
                            "RG.dat")




@dataclass(frozen=True)
class CaptureResult:
	test_dir: Path
	source_file: Path
	final_line: str


@dataclass(frozen=True)
class WriteResult:
	source_filename: str
	output_file: Path
	row_count: int


def _extract_test_number(test_dir_name: str) -> int | None:
	"""Parse a simulation test directory name into its numeric index.

	The script expects test case folders to follow the naming pattern
	``test_<number>`` (for example ``test_1`` or ``test_42``). This helper
	validates that format and extracts the integer component for sorting.

	Parameters
	----------
	test_dir_name : str
		Directory name to evaluate.

	Returns
	-------
	int | None
		Parsed integer test number when the name matches ``test_<number>``;
		otherwise ``None``.
	"""
	match = re.fullmatch(r"test_(\d+)", test_dir_name)
	if not match:
		return None
	return int(match.group(1))


def _discover_test_dirs(tests_root: Path) -> list[Path]:
	"""Find and numerically sort valid ``test_*`` directories.

	Only immediate child directories under ``tests_root`` are considered.
	Entries that do not match the ``test_<number>`` naming scheme are ignored.

	Parameters
	----------
	tests_root : Path
		Root directory containing simulation test folders.

	Returns
	-------
	list[Path]
		List of matching test directories sorted by numeric test index.
	"""
	candidates = []
	for path in tests_root.iterdir():
		if not path.is_dir():
			continue
		test_num = _extract_test_number(path.name)
		if test_num is None:
			continue
		candidates.append((test_num, path))

	return [path for _, path in sorted(candidates, key=lambda pair: pair[0])]


def _read_final_nonempty_line(path: Path) -> str:
	"""Return the final non-blank line from a text file.

	The file is read in full, each line is stripped, and empty lines are
	ignored. The last remaining line is used as the expected output value.

	Parameters
	----------
	path : Path
		Text file to inspect.

	Returns
	-------
	str
		Last non-empty, stripped line in the file.

	Raises
	------
	ValueError
		Raised when the file contains no non-empty lines.
	"""
	lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
	if not lines:
		raise ValueError(f"File has no non-empty lines: {path}")
	return lines[-1]


def _combined_output_name(source_filename: str) -> str:
	"""Build the output filename for a combined expected-output artifact.

	Nested path separators are replaced to keep the result filesystem-safe,
	then the ``.final_lines.txt`` suffix is appended.

	Parameters
	----------
	source_filename : str
		Original filename captured from each test directory.

	Returns
	-------
	str
		Filename used for the aggregated expected-output file.
	"""
	# Keep names filesystem-safe if nested paths are ever passed via --files.
	safe_name = source_filename.replace("/", "__")
	return f"{safe_name}.final_lines.txt"


def _write_combined_outputs(
	expected_root: Path,
	results: list[CaptureResult],
) -> list[WriteResult]:
	"""Write one consolidated output file per source filename.

	Input capture results are grouped by the original source file name
	(e.g. ``ENERGY.dat``). For each group, rows are written as tab-separated
	``test_dir_name`` and ``final_line`` values in test-number order.

	Parameters
	----------
	expected_root : Path
		Directory where consolidated expected-output files are written.
	results : list[CaptureResult]
		Collected final-line values for each test/file pair.

	Returns
	-------
	list[WriteResult]
		Per-file write metadata including destination path and row count.
	"""
	grouped: dict[str, list[CaptureResult]] = {}
	for result in results:
		filename = result.source_file.name
		grouped.setdefault(filename, []).append(result)

	write_results: list[WriteResult] = []
	for source_filename in sorted(grouped):
		group_results = sorted(
			grouped[source_filename],
			key=lambda item: _extract_test_number(item.test_dir.name) or 0,
		)
		output_path = expected_root / _combined_output_name(source_filename)
		lines = [f"{item.test_dir.name}\t{item.final_line}" for item in group_results]
		output_path.write_text("\n".join(lines) + "\n")

		write_results.append(
			WriteResult(
				source_filename=source_filename,
				output_file=output_path,
				row_count=len(group_results),
			)
		)

	return write_results


def _capture_expected_outputs(
	tests_root: Path,
	expected_root: Path,
	files_to_capture: list[str],
	allow_missing: bool,
) -> tuple[list[CaptureResult], list[WriteResult], list[str]]:
	"""Collect final lines from simulation outputs and write grouped artifacts.

	For each discovered ``test_*`` directory, this function attempts to read
	the final non-empty line from each requested source file. Missing files
	and empty-file conditions are recorded as issues. Successful captures are
	then aggregated into per-source combined output files.

	Parameters
	----------
	tests_root : Path
		Directory containing ``test_*`` folders.
	expected_root : Path
		Directory where combined expected-output files will be written.
	files_to_capture : list[str]
		Filenames to process inside each test directory.
	allow_missing : bool
		If ``True``, missing files are emitted as warnings; if ``False``,
		they are emitted as errors.

	Returns
	-------
	tuple[list[CaptureResult], list[WriteResult], list[str]]
		Tuple of:
		1. Successful capture records.
		2. Metadata for combined output files that were written.
		3. Issue messages prefixed with ``[WARN]`` or ``[ERROR]``.
	"""
	expected_root.mkdir(parents=True, exist_ok=True)

	results: list[CaptureResult] = []
	write_results: list[WriteResult] = []
	issues: list[str] = []

	test_dirs = _discover_test_dirs(tests_root)
	if not test_dirs:
		issues.append(f"No test directories found in {tests_root}")
		return results, write_results, issues

	for test_dir in test_dirs:
		for filename in files_to_capture:
			source_path = test_dir / filename
			if not source_path.exists():
				message = f"Missing file for {test_dir.name}: {filename}"
				if allow_missing:
					issues.append(f"[WARN] {message}")
					continue
				issues.append(f"[ERROR] {message}")
				continue

			try:
				final_line = _read_final_nonempty_line(source_path)
			except ValueError as exc:
				issues.append(f"[ERROR] {exc}")
				continue

			results.append(
				CaptureResult(
					test_dir=test_dir,
					source_file=source_path,
					final_line=final_line,
				)
			)

	write_results = _write_combined_outputs(expected_root=expected_root, results=results)
	return results, write_results, issues


def _build_arg_parser() -> argparse.ArgumentParser:
	"""Create the command-line argument parser for this utility.

	The parser exposes options for selecting test directories, choosing the
	output destination, overriding files to inspect, and controlling whether
	missing files should fail the command.

	Returns
	-------
	argparse.ArgumentParser
		Configured parser for the expected-output generation CLI.
	"""
	parser = argparse.ArgumentParser(
		description=(
			"Generate stand-alone expected-output files by taking the final non-empty "
			"line from selected simulation output files in each test_* directory."
		)
	)
	parser.add_argument(
		"--tests-root",
		type=Path,
		default=Path(__file__).resolve().parent,
		help="Directory containing test_1, test_2, ... folders (default: current test suite folder).",
	)
	parser.add_argument(
		"--expected-root",
		type=Path,
		default=None,
		help="Output directory for expected-output files (default: <tests-root>/expected_output).",
	)
	parser.add_argument(
		"--files",
		nargs="+",
		default=list(DEFAULT_FILES_TO_CAPTURE),
		help=(
			"List of file names to inspect inside each test directory. "
			"Default: {}".format(" ".join(DEFAULT_FILES_TO_CAPTURE))
		),
	)
	parser.add_argument(
		"--allow-missing",
		action="store_true",
		help="Skip missing files with warnings instead of returning an error code.",
	)
	return parser


def main(argv: list[str] | None = None) -> int:
	"""Run the expected-output generation CLI workflow.

	This entry point parses command-line arguments, resolves paths, captures
	final lines from requested simulation output files, writes combined files,
	prints a summary, and returns a process exit code.

	Parameters
	----------
	argv : list[str] | None, optional
		Argument vector to parse. When ``None``, arguments are read from
		``sys.argv`` by ``argparse``.

	Returns
	-------
	int
		Exit status code: ``0`` on success, ``1`` if any ``[ERROR]`` issues
		were recorded.
	"""
	parser = _build_arg_parser()
	args = parser.parse_args(argv)

	tests_root = args.tests_root.resolve()
	expected_root = (
		args.expected_root.resolve()
		if args.expected_root is not None
		else (tests_root / "expected_output").resolve()
	)

	results, write_results, issues = _capture_expected_outputs(
		tests_root=tests_root,
		expected_root=expected_root,
		files_to_capture=args.files,
		allow_missing=args.allow_missing,
	)

	for write_result in write_results:
		print(
			f"[OK] {write_result.source_filename} -> {write_result.output_file.name} "
			f"({write_result.row_count} tests)"
		)

	for issue in issues:
		print(issue)

	has_error = any(issue.startswith("[ERROR]") for issue in issues)
	print(
		f"Captured {len(results)} final lines and wrote {len(write_results)} "
		f"combined expected-output files to {expected_root}"
	)
	return 1 if has_error else 0


if __name__ == "__main__":
	sys.exit(main())
