#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np


DEFAULT_DATA_FILES = [
    "eedf",
    "swarmParameters",
    "rateCoefficients",
    "powerBalance",
    "lookUpTable",
]


@dataclass(frozen=True)
class RunSpec:
    cells: int
    is_reference: bool


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run a LoKI-B convergence study by sweeping the energy-grid cell count, "
            "writing one labeled JSON result file per run."
        )
    )
    parser.add_argument("input_file", type=Path, help="Path to the LoKI input file (.json or .in).")
    parser.add_argument("min_cells", type=int, help="Minimum number of cells to simulate.")
    parser.add_argument("max_cells", type=int, help="Maximum number of cells to simulate.")
    parser.add_argument("step_size", type=int, help="Cell-count increment between runs.")
    parser.add_argument(
        "output_dir",
        type=Path,
        help="Directory where the labeled JSON result files and the manifest will be stored.",
    )
    parser.add_argument(
        "--reference-cells",
        type=int,
        default=None,
        help=(
            "Optional extra high-resolution reference run. If omitted, the largest run in the study "
            "is used as the reference."
        ),
    )
    parser.add_argument(
        "--loki-bin",
        type=Path,
        default=None,
        help="Path to the loki executable. Defaults to build/app/loki or loki from PATH.",
    )
    parser.add_argument(
        "--legacytojson-bin",
        type=Path,
        default=None,
        help=(
            "Path to the loki_legacytojson executable. Defaults to build/app/loki_legacytojson "
            "or loki_legacytojson from PATH."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite an existing manifest or result files in the output directory.",
    )
    return parser.parse_args()


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def resolve_binary(explicit_path: Path | None, build_name: str, fallback_name: str, option_name: str) -> Path:
    if explicit_path is not None:
        path = explicit_path.expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Executable not found: {path}")
        return path

    build_path = repo_root() / "build" / "app" / build_name
    if build_path.is_file():
        return build_path

    discovered = shutil.which(fallback_name)
    if discovered is not None:
        return Path(discovered).resolve()

    raise FileNotFoundError(
        f"Could not find '{fallback_name}'. Pass {option_name} explicitly or build LoKI first."
    )


def run_command(command: list[str], *, cwd: Path, capture_stdout: bool) -> str:
    result = subprocess.run(
        command,
        cwd=cwd,
        stdout=subprocess.PIPE if capture_stdout else subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        message = [
            f"Command failed with exit code {result.returncode}: {' '.join(command)}",
            "--- stdout ---",
            result.stdout.rstrip(),
            "--- stderr ---",
            result.stderr.rstrip(),
        ]
        raise RuntimeError("\n".join(message))
    return result.stdout if result.stdout is not None else ""


def needs_legacy_json_conversion(config: dict[str, Any]) -> bool:
    working_conditions = config.get("workingConditions")
    if isinstance(working_conditions, dict):
        for key in (
            "gasPressure",
            "gasTemperature",
            "electronDensity",
            "electronTemperature",
            "excitationFrequency",
            "reducedField",
        ):
            value = working_conditions.get(key)
            if isinstance(value, (int, float, str)):
                return True

    electron_kinetics = config.get("electronKinetics")
    if not isinstance(electron_kinetics, dict):
        return False

    gas_properties = electron_kinetics.get("gasProperties")
    if isinstance(gas_properties, dict) and isinstance(gas_properties.get("fraction"), list):
        return True

    state_properties = electron_kinetics.get("stateProperties")
    if isinstance(state_properties, dict):
        for value in state_properties.values():
            if isinstance(value, list):
                return True

    return False


def load_config(input_path: Path, legacytojson_bin: Path) -> dict[str, Any]:
    if input_path.suffix == ".in":
        converted = run_command([str(legacytojson_bin), str(input_path)], cwd=input_path.parent, capture_stdout=True)
        return json.loads(converted)

    if input_path.suffix != ".json":
        raise ValueError(f"Unsupported input format '{input_path.suffix}'. Expected .json or .in.")

    with input_path.open(encoding="utf-8") as handle:
        config = json.load(handle)

    if needs_legacy_json_conversion(config):
        converted = run_command([str(legacytojson_bin), str(input_path)], cwd=input_path.parent, capture_stdout=True)
        return json.loads(converted)

    return config


def ensure_output_section(config: dict[str, Any], json_output_path: Path) -> None:
    output = config.setdefault("output", {})
    if not isinstance(output, dict):
        raise ValueError("Expected 'output' to be an object in the input configuration.")

    output["isOn"] = True
    output["writeJSON"] = True
    output["writeText"] = False
    output["JSONFile"] = str(json_output_path)
    output.setdefault("folder", json_output_path.stem)

    data_files = output.get("dataFiles")
    if data_files is None:
        output["dataFiles"] = DEFAULT_DATA_FILES.copy()
        return

    if not isinstance(data_files, list):
        raise ValueError("Expected 'output.dataFiles' to be an array when present.")

    merged = list(dict.fromkeys([*data_files, "swarmParameters"]))
    output["dataFiles"] = merged


def set_cell_count(config: dict[str, Any], cells: int) -> None:
    try:
        energy_grid = config["electronKinetics"]["numerics"]["energyGrid"]
    except KeyError as exc:
        raise KeyError("Could not find electronKinetics.numerics.energyGrid in the input configuration.") from exc

    if not isinstance(energy_grid, dict):
        raise ValueError("Expected electronKinetics.numerics.energyGrid to be an object.")

    config["electronKinetics"]["numerics"]["energyGrid"] = {
        "nonuniformGrid": {
            "maxEnergy": energy_grid["maxEnergy"],
            "nodeDistribution": list(np.linspace(0.0, 1.0, cells + 1) ** 1.5)
        }
    }
    return

    # if "cellNumber" in energy_grid:
    #     energy_grid["cellNumber"] = cells
    #     return

    # nonuniform_grid = energy_grid.get("nonuniformGrid")
    # if isinstance(nonuniform_grid, dict):
    #     node_distribution = nonuniform_grid.get("nodeDistribution")
    #     if isinstance(node_distribution, dict) and "nvalues" in node_distribution:
    #         node_distribution["nvalues"] = cells + 1
    #         return
    #     raise ValueError(
    #         "This convergence helper only supports nonuniform grids whose nodeDistribution defines 'nvalues'."
    #     )

    raise ValueError("Could not map the requested cell count onto the energy-grid configuration.")


def write_temp_config(input_path: Path, config: dict[str, Any], cells: int) -> Path:
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        dir=input_path.parent,
        prefix=f".{input_path.stem}_cells_{cells}_",
        suffix=".json",
        delete=False,
    ) as handle:
        json.dump(config, handle, indent=2)
        handle.write("\n")
        return Path(handle.name)


def build_run_specs(min_cells: int, max_cells: int, step_size: int, reference_cells: int | None) -> list[RunSpec]:
    if min_cells <= 0 or max_cells <= 0 or step_size <= 0:
        raise ValueError("min_cells, max_cells, and step_size must all be positive integers.")
    if min_cells > max_cells:
        raise ValueError("min_cells must not be larger than max_cells.")
    if reference_cells is not None and reference_cells <= 0:
        raise ValueError("reference_cells must be positive when provided.")

    study_cells = list(range(min_cells, max_cells + 1, step_size))
    if not study_cells:
        raise ValueError("The requested cell range did not generate any runs.")

    effective_reference = reference_cells if reference_cells is not None else study_cells[-1]

    runs = []
    for cells in sorted(set([*study_cells, effective_reference])):
        runs.append(RunSpec(cells=cells, is_reference=(cells == effective_reference)))
    return runs


def annotate_result_file(
    result_path: Path,
    *,
    input_path: Path,
    cells: int,
    reference_cells: int,
    is_reference: bool,
) -> None:
    with result_path.open(encoding="utf-8") as handle:
        data = json.load(handle)

    data["convergenceStudy"] = {
        "sourceInputFile": str(input_path.resolve()),
        "cellNumber": cells,
        "referenceCellNumber": reference_cells,
        "isReference": is_reference,
    }

    with result_path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2)
        handle.write("\n")


def main() -> int:
    args = parse_args()

    input_path = args.input_file.expanduser().resolve()
    if not input_path.is_file():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    loki_bin = resolve_binary(args.loki_bin, "loki", "loki", "--loki-bin")
    legacytojson_bin = resolve_binary(
        args.legacytojson_bin,
        "loki_legacytojson",
        "loki_legacytojson",
        "--legacytojson-bin",
    )

    run_specs = build_run_specs(args.min_cells, args.max_cells, args.step_size, args.reference_cells)
    reference_cells = next(spec.cells for spec in run_specs if spec.is_reference)
    study_cells = list(range(args.min_cells, args.max_cells + 1, args.step_size))

    manifest_path = output_dir / "convergence_manifest.json"
    if manifest_path.exists() and not args.force:
        raise FileExistsError(f"Manifest already exists: {manifest_path}. Re-run with --force to overwrite.")

    width = max(4, max(len(str(spec.cells)) for spec in run_specs))
    base_config = load_config(input_path, legacytojson_bin)

    manifest: dict[str, Any] = {
        "sourceInputFile": str(input_path),
        "lokiBinary": str(loki_bin),
        "legacyToJsonBinary": str(legacytojson_bin),
        "studyCellNumbers": study_cells,
        "referenceCellNumber": reference_cells,
        "runs": [],
    }

    for spec in run_specs:
        result_path = output_dir / f"cells_{spec.cells:0{width}d}.json"
        if result_path.exists() and not args.force:
            raise FileExistsError(f"Result file already exists: {result_path}. Re-run with --force to overwrite.")

        config = json.loads(json.dumps(base_config))
        set_cell_count(config, spec.cells)
        ensure_output_section(config, result_path)
        temp_config = write_temp_config(input_path, config, spec.cells)

        try:
            run_command([str(loki_bin), str(temp_config)], cwd=input_path.parent, capture_stdout=False)
        finally:
            temp_config.unlink(missing_ok=True)

        annotate_result_file(
            result_path,
            input_path=input_path,
            cells=spec.cells,
            reference_cells=reference_cells,
            is_reference=spec.is_reference,
        )

        manifest["runs"].append(
            {
                "cellNumber": spec.cells,
                "jsonFile": result_path.name,
                "isReference": spec.is_reference,
            }
        )

    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2)
        handle.write("\n")

    print(f"Wrote {len(run_specs)} JSON result file(s) to {output_dir}")
    print(f"Reference run: {reference_cells} cells")
    print(f"Manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        raise SystemExit(1)
