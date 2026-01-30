from __future__ import annotations

from pathlib import Path
import argparse
import configparser


def _num(line: str) -> float:
    for token in line.replace(",", " ").split():
        try:
            return float(token)
        except ValueError:
            continue
    raise ValueError(f"Expected numeric value in line: {line!r}")


def parse_def_file(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if len(lines) < 34:
        raise ValueError(f"Definition file too short: {path}")

    idx = 0
    _ = lines[idx]  # description
    idx += 1
    input_file = lines[idx].strip().split()[0]
    idx += 1
    output_file = lines[idx].strip().split()[0]
    idx += 1

    nums = [_num(line) for line in lines[idx:] if line.strip()]

    if len(nums) == 33:
        hour_start_input, hour_start_output, hour_end_output = nums[:3]
        cursor = 3
    elif len(nums) == 31:
        hour_start_input = nums[0]
        hour_start_output = nums[0]
        hour_end_output = nums[0]
        cursor = 1
    else:
        raise ValueError(f"Unexpected numeric line count ({len(nums)}) in {path}")

    hour_start = hour_start_output
    hz = nums[cursor]; cursor += 1
    hz_out = nums[cursor]; cursor += 1
    pr_offset = nums[cursor]; cursor += 1
    salinity = nums[cursor]; cursor += 1
    flag_pr = int(nums[cursor]); cursor += 1
    elevation_m = nums[cursor]; cursor += 1
    atmpr_mbar = nums[cursor]; cursor += 1
    at = nums[cursor]; cursor += 1
    bt = nums[cursor]; cursor += 1
    ct = nums[cursor]; cursor += 1
    dt = nums[cursor]; cursor += 1
    ao2 = nums[cursor]; cursor += 1
    bo2 = nums[cursor]; cursor += 1
    co2 = nums[cursor]; cursor += 1
    do2 = nums[cursor]; cursor += 1
    eo2 = nums[cursor]; cursor += 1  # kept for completeness
    fo2 = nums[cursor]; cursor += 1
    go2 = nums[cursor]; cursor += 1
    ho2 = nums[cursor]; cursor += 1
    _flahout = nums[cursor]; cursor += 1  # unused in Python port
    time1 = nums[cursor]; cursor += 1
    time2 = nums[cursor]; cursor += 1
    time3 = nums[cursor]; cursor += 1
    o2factor1 = nums[cursor]; cursor += 1
    o2factor2 = nums[cursor]; cursor += 1
    o2factor3 = nums[cursor]; cursor += 1
    o24 = nums[cursor]; cursor += 1
    o25 = nums[cursor]; cursor += 1
    o2factor4 = nums[cursor]; cursor += 1
    o2factor5 = nums[cursor]; cursor += 1

    return {
        "input_file": input_file,
        "output_file": output_file,
        "hour_start_input": hour_start_input,
        "hour_start_output": hour_start_output,
        "hour_end_output": hour_end_output,
        "hour_start": hour_start,
        "hz": hz,
        "hz_out": hz_out,
        "pr_offset": pr_offset,
        "salinity": salinity,
        "flag_pr": flag_pr,
        "elevation_m": elevation_m,
        "atmpr_mbar": atmpr_mbar,
        "temp_coeffs": [at, bt, ct, dt],
        "o2_coeffs": [ao2, bo2, co2, do2, eo2, fo2, go2, ho2],
        "time_points": [time1, time2, time3],
        "time_factors": [o2factor1, o2factor2, o2factor3],
        "o2_points": [o24, o25],
        "o2_factors": [o2factor4, o2factor5],
    }


def write_ini(data: dict, output_path: Path) -> None:
    config = configparser.ConfigParser()
    config["paths"] = {
        "input_file": str(data["input_file"]),
        "output_file": str(data["output_file"]),
    }
    config["sampling"] = {
        "hour_start_input": str(data["hour_start_input"]),
        "hour_start_output": str(data["hour_start_output"]),
        "hour_end_output": str(data["hour_end_output"]),
        "hour_start": str(data["hour_start"]),
        "hz": str(data["hz"]),
        "hz_out": str(data["hz_out"]),
    }
    config["environment"] = {
        "pr_offset": str(data["pr_offset"]),
        "salinity": str(data["salinity"]),
        "flag_pr": str(data["flag_pr"]),
        "elevation_m": str(data["elevation_m"]),
        "atmpr_mbar": str(data["atmpr_mbar"]),
    }
    config["temp_coeffs"] = {
        "coeffs": ", ".join(str(v) for v in data["temp_coeffs"]),
    }
    config["o2_coeffs"] = {
        "coeffs": ", ".join(str(v) for v in data["o2_coeffs"]),
    }
    config["calibration_time"] = {
        "time_points": ", ".join(str(v) for v in data["time_points"]),
        "factors": ", ".join(str(v) for v in data["time_factors"]),
    }
    config["calibration_concentration"] = {
        "o2_points": ", ".join(str(v) for v in data["o2_points"]),
        "factors": ", ".join(str(v) for v in data["o2_factors"]),
    }

    with output_path.open("w", encoding="utf-8") as handle:
        config.write(handle)


def main() -> int:
    parser = argparse.ArgumentParser(description="Convert ExtractRINKO .def to vector_formatter.ini")
    parser.add_argument("def_path", type=Path, help="Path to the .def file")
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output .ini path (default: same path with .ini extension)",
    )
    args = parser.parse_args()

    output_path = args.out or args.def_path.with_suffix(".ini")
    data = parse_def_file(args.def_path)
    write_ini(data, output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
