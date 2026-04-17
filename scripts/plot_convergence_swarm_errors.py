#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median
from typing import Any


PARAMETERS_TO_SKIP = {"Reduced electric field"}
RESULT_FILE_RE = re.compile(r"cells_(\d+)\.json$")


@dataclass(frozen=True)
class RunData:
    cells: int
    path: Path
    jobs: dict[str, dict[str, float]]
    is_reference: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot the absolute relative error of LoKI-B swarm parameters from a convergence study."
        )
    )
    parser.add_argument(
        "study_dir",
        type=Path,
        help="Directory containing convergence JSON files, typically created by run_convergence_study.py.",
    )
    parser.add_argument(
        "--reference-file",
        type=Path,
        default=None,
        help="Use a specific JSON result file as the reference instead of the default manifest or largest-cell run.",
    )
    parser.add_argument(
        "--reference-cells",
        type=int,
        default=None,
        help="Use the run with this many cells as the reference.",
    )
    parser.add_argument(
        "--aggregate",
        choices=("mean", "max", "median"),
        default="max",
        help="How to reduce job-by-job errors into one error per parameter and cell count.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Path of the output figure. Defaults to <study_dir>/swarm_parameter_relative_error.png.",
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        default=None,
        help="Path of the detailed error summary JSON. Defaults to <study_dir>/swarm_parameter_relative_error.json.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure interactively after writing it.",
    )
    return parser.parse_args()


def load_json(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return json.load(handle)


def parse_cells_from_name(path: Path) -> int:
    match = RESULT_FILE_RE.fullmatch(path.name)
    if match is None:
        raise ValueError(f"Could not infer the cell count from file name '{path.name}'.")
    return int(match.group(1))


def parse_swarm_parameters(entries: Any) -> dict[str, float]:
    if not isinstance(entries, list):
        raise ValueError("Expected 'swarm_parameters' to be an array.")

    result: dict[str, float] = {}
    for entry in entries:
        if not isinstance(entry, dict) or len(entry) != 1:
            raise ValueError("Unexpected swarm parameter entry structure.")
        name, payload = next(iter(entry.items()))
        if not isinstance(payload, dict) or "value" not in payload:
            raise ValueError(f"Missing value for swarm parameter '{name}'.")
        result[name] = float(payload["value"])
    return result


def extract_jobs(data: dict[str, Any]) -> dict[str, dict[str, float]]:
    jobs: dict[str, dict[str, float]] = {}
    for key, value in data.items():
        if not isinstance(value, dict):
            continue
        if "swarm_parameters" not in value:
            continue
        jobs[key] = parse_swarm_parameters(value["swarm_parameters"])

    if not jobs:
        raise ValueError("No top-level simulation entries with 'swarm_parameters' were found.")
    return jobs


def load_run(path: Path) -> RunData:
    data = load_json(path)
    metadata = data.get("convergenceStudy", {})
    cells = int(metadata.get("cellNumber", parse_cells_from_name(path)))
    is_reference = bool(metadata.get("isReference", False))
    return RunData(cells=cells, path=path, jobs=extract_jobs(data), is_reference=is_reference)


def discover_runs(study_dir: Path) -> list[RunData]:
    manifest_path = study_dir / "convergence_manifest.json"
    if manifest_path.is_file():
        manifest = load_json(manifest_path)
        runs = [
            load_run(study_dir / run["jsonFile"])
            for run in manifest.get("runs", [])
        ]
        if runs:
            return sorted(runs, key=lambda run: run.cells)

    runs = [load_run(path) for path in sorted(study_dir.glob("cells_*.json"))]
    if not runs:
        raise FileNotFoundError(f"No convergence result files were found in {study_dir}.")
    return sorted(runs, key=lambda run: run.cells)


def choose_reference(
    runs: list[RunData],
    reference_file: Path | None,
    reference_cells: int | None,
) -> RunData:
    if reference_file is not None:
        target = reference_file.expanduser().resolve()
        for run in runs:
            if run.path.resolve() == target:
                return run
        return load_run(target)

    if reference_cells is not None:
        for run in runs:
            if run.cells == reference_cells:
                return run
        raise ValueError(f"No run with {reference_cells} cells was found.")

    marked = [run for run in runs if run.is_reference]
    if marked:
        return marked[0]

    return max(runs, key=lambda run: run.cells)


def relative_error(value: float, reference: float) -> float:
    if math.isclose(reference, 0.0, abs_tol=1e-30):
        if math.isclose(value, 0.0, abs_tol=1e-30):
            return 0.0
        return math.inf
    return abs((value - reference) / reference)


def aggregate(values: list[float], mode: str) -> float:
    if any(math.isinf(value) for value in values):
        return math.inf
    finite_values = [value for value in values if math.isfinite(value)]
    if not finite_values:
        return math.inf
    if mode == "mean":
        return mean(finite_values)
    if mode == "median":
        return median(finite_values)
    if mode == "max":
        return max(finite_values)
    raise ValueError(f"Unsupported aggregation mode '{mode}'.")


def compute_errors(runs: list[RunData], reference: RunData, mode: str) -> dict[str, Any]:
    common_jobs = set(reference.jobs)
    for run in runs:
        common_jobs &= set(run.jobs)
    if not common_jobs:
        raise ValueError("The selected runs do not share any common simulation jobs.")

    ordered_jobs = sorted(common_jobs)
    parameter_names = [
        name
        for name in reference.jobs[ordered_jobs[0]]
        if name not in PARAMETERS_TO_SKIP
    ]

    summary_runs = []
    aggregate_by_parameter = {name: [] for name in parameter_names}

    for run in sorted(runs, key=lambda item: item.cells):
        if run.is_reference == True:
            continue

        run_entry: dict[str, Any] = {
            "cellNumber": run.cells,
            "jsonFile": run.path.name,
            "isReference": run.path.resolve() == reference.path.resolve(),
            "aggregateErrors": {},
            "jobErrors": {},
        }

        for parameter in parameter_names:
            job_errors: dict[str, float] = {}
            values: list[float] = []
            for job in ordered_jobs:
                error = relative_error(run.jobs[job][parameter], reference.jobs[job][parameter])
                job_errors[job] = error
                values.append(error)

            aggregated = aggregate(values, mode)
            run_entry["jobErrors"][parameter] = job_errors
            run_entry["aggregateErrors"][parameter] = aggregated
            aggregate_by_parameter[parameter].append((run.cells, aggregated))

        summary_runs.append(run_entry)

    return {
        "referenceFile": reference.path.name,
        "referenceCellNumber": reference.cells,
        "aggregateMode": mode,
        "commonJobs": ordered_jobs,
        "parameters": parameter_names,
        "runs": summary_runs,
        "aggregateByParameter": {
            parameter: [
                {"cellNumber": cells, "error": error}
                for cells, error in values
            ]
            for parameter, values in aggregate_by_parameter.items()
        },
    }


def import_matplotlib() -> Any:
    try:
        import matplotlib.pyplot as plt
    except ImportError as exc:
        raise RuntimeError("matplotlib is required for plotting. Install it with 'pip install matplotlib'.") from exc
    return plt


def plot_summary(summary: dict[str, Any], output_path: Path, show: bool) -> None:
    plt = import_matplotlib()

    fig, ax = plt.subplots(figsize=(10, 6))
    has_finite_data = False

    for parameter, series in summary["aggregateByParameter"].items():
        x_values = [point["cellNumber"] for point in series]
        y_values = [point["error"] for point in series]
        finite_y = [value for value in y_values if math.isfinite(value) and value > 0.0]
        has_finite_data = has_finite_data or any(math.isfinite(value) for value in y_values)

        plot_y = [math.nan if not math.isfinite(value) else value for value in y_values]
        ax.plot(x_values, plot_y, marker="o", label=parameter)

        if finite_y:
            ax.set_yscale("log")

    ax.set_xlabel("Number of cells")
    ax.set_ylabel(f"Absolute relative error ({summary['aggregateMode']})")
    ax.set_title(
        f"Swarm-parameter convergence relative to {summary['referenceCellNumber']} cells"
    )
    ax.grid(True, which="both", linestyle="--", linewidth=0.5, alpha=0.7)
    ax.legend()
    fig.tight_layout()

    if not has_finite_data:
        raise RuntimeError("No finite errors were available to plot.")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    if show:
        plt.show()
    plt.close(fig)


def main() -> int:
    args = parse_args()

    study_dir = args.study_dir.expanduser().resolve()
    if not study_dir.is_dir():
        raise FileNotFoundError(f"Study directory not found: {study_dir}")

    runs = discover_runs(study_dir)
    reference = choose_reference(runs, args.reference_file, args.reference_cells)
    summary = compute_errors(runs, reference, args.aggregate)

    output_path = (
        args.output.expanduser().resolve()
        if args.output is not None
        else study_dir / "swarm_parameter_relative_error.png"
    )
    summary_json_path = (
        args.summary_json.expanduser().resolve()
        if args.summary_json is not None
        else study_dir / "swarm_parameter_relative_error.json"
    )

    with summary_json_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    plot_summary(summary, output_path, args.show)

    print(f"Wrote plot to {output_path}")
    print(f"Wrote detailed error summary to {summary_json_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1)
