from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import argparse
import configparser
import math
import re
try:
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover - lightweight fallback when numpy isn't available
    class _Poly1D:
        def __init__(self, coeffs: Sequence[float]) -> None:
            self.coeffs = list(coeffs)

        def __call__(self, x: float) -> float:
            acc = 0.0
            for coeff in self.coeffs:
                acc = acc * x + coeff
            return acc

    class _NumpyFallback:
        @staticmethod
        def poly1d(coeffs: Sequence[float]) -> _Poly1D:
            return _Poly1D(coeffs)

    np = _NumpyFallback()
from typing import Dict, Iterable, List, Sequence


@dataclass(frozen=True)
class ColumnMap:
    # 1-based column indexes in the converted Nortek Vector text file.
    vx: int = 3
    vy: int = 4
    vz: int = 5
    pr: int = 15
    p_counts: int = 16  # oxygen counts
    q_counts: int = 17  # temperature counts


@dataclass(frozen=True)
class Definition:
    input_file: str
    output_file: str
    hour_start_input: float
    hour_start: float
    hour_end_output: float
    hz: float
    hz_out: float
    pr_offset: float
    salinity: float
    flag_pr: int
    elevation_m: float
    atmpr_mbar: float
    temp_coeffs: Sequence[float]  # at, bt, ct, dt
    o2_coeffs: Sequence[float]  # ao2, bo2, co2, do2, eo2, fo2, go2, ho2
    time_points: Sequence[float]  # time1, time2, time3
    time_factors: Sequence[float]  # o2factor1, o2factor2, o2factor3
    o2_points: Sequence[float]  # o24, o25
    o2_factors: Sequence[float]  # o2factor4, o2factor5


# Hard-coded output columns (edit when the final column set is defined).
OUTPUT_COLUMNS: Sequence[str] = (
    "hour",
    "vx",
    "vy",
    "vz",
    "pr_mbar",
    "temp_c",
    "o2_umol_l",
    "o2sat_umol_l",
)


DEFAULT_FLOAT_FORMAT = "{:.6f}"


def _float_list(raw: str, *, expected: int | None = None) -> List[float]:
    tokens = [tok for tok in re.split(r"[,\s]+", raw.strip()) if tok]
    values = [float(tok) for tok in tokens]
    if expected is not None and len(values) != expected:
        raise ValueError(f"Expected {expected} values, got {len(values)}: {raw!r}")
    return values


def parse_ini_file(path: Path) -> Definition:
    # Parse the vector_formatter.ini configuration.
    config = configparser.ConfigParser()
    config.read(path, encoding="utf-8")

    def _section(name: str) -> configparser.SectionProxy:
        if name not in config:
            raise ValueError(f"Missing [{name}] section in {path}")
        return config[name]

    paths = _section("paths")
    sampling = _section("sampling")
    environment = _section("environment")
    temp_coeffs_section = _section("temp_coeffs")
    o2_coeffs_section = _section("o2_coeffs")
    time_cal = _section("calibration_time")
    conc_cal = _section("calibration_concentration")

    input_file = paths.get("input_file", "")
    output_file = paths.get("output_file", "").strip()
    if not output_file:
        default_name = f"{Path(input_file).stem}_formatted.dat"
        output_file = default_name

    hour_start = sampling.getfloat("hour_start_output", fallback=sampling.getfloat("hour_start"))
    hour_start_input = sampling.getfloat("hour_start_input", fallback=hour_start)
    hour_end_output = sampling.getfloat("hour_end_output", fallback=hour_start)
    hz = sampling.getfloat("hz")
    hz_out = sampling.getfloat("hz_out")

    pr_offset = environment.getfloat("pr_offset")
    salinity = environment.getfloat("salinity")
    flag_pr = environment.getint("flag_pr")
    elevation_m = environment.getfloat("elevation_m")
    atmpr_mbar = environment.getfloat("atmpr_mbar")

    temp_coeffs = _float_list(temp_coeffs_section.get("coeffs", ""), expected=4)
    o2_coeffs = _float_list(o2_coeffs_section.get("coeffs", ""), expected=8)

    time_points = _float_list(time_cal.get("time_points", ""), expected=3)
    time_factors = _float_list(time_cal.get("factors", ""), expected=3)
    o2_points = _float_list(conc_cal.get("o2_points", ""), expected=2)
    o2_factors = _float_list(conc_cal.get("factors", ""), expected=2)

    fields = Definition.__dataclass_fields__.keys()
    return Definition(**{name: locals()[name] for name in fields})


def o2sat(temp_c: float, sal: float) -> float:
    # Garcia and Gordon (1992), oxygen saturation in umol kg-1.
    a_coeffs = [1.41567, 1.01567, 4.93845, 4.11890, 3.20684, 5.80818]  # a5..a0
    b_coeffs = [-5.54491e-03, -7.93334e-03, -7.25958e-03, -7.01211e-03]  # b3..b0
    c0 = -1.32412e-07

    ts = math.log((298.15 - temp_c) / (273.15 + temp_c))
    poly_a = np.poly1d(a_coeffs)(ts)
    poly_b = np.poly1d(b_coeffs)(ts)
    return math.exp(poly_a + sal * poly_b + c0 * sal ** 2)


def density(temp_c: float, sal: float) -> float:
    # Standard Methods for the Examination of Water and Wastewater, density in g cm-3.
    a_coeffs = [5.3875e-09, -8.2467e-07, 7.6438e-05, -4.0899e-03, 8.24493e-01]
    b_coeffs = [-1.6546e-06, 1.0227e-04, -5.72466e-03]
    c = 4.83140e-04
    density0_coeffs = [6.536332e-09, -1.120083e-06, 1.001685e-04, -9.095290e-03, 6.793952e-02, 999.842594]
    a = np.poly1d(a_coeffs)(temp_c)
    b = np.poly1d(b_coeffs)(temp_c)
    density0 = np.poly1d(density0_coeffs)(temp_c)
    return (density0 + a * sal + b * sal ** 1.5 + c * sal ** 2) / 1000.0


def o2_rinko(o2_coeffs: Sequence[float], temp_c: float, volt: float, o2sat_den: float) -> float:
    # Compute oxygen concentration (umol L-1) from the RINKO sensor model.
    ao2, bo2, co2, do2, _eo2, fo2, go2, ho2 = o2_coeffs
    dt = temp_c - 25.0
    denom = np.poly1d([fo2, do2, 1.0])(dt)
    o2pct = go2 + ho2 * (ao2 / denom + bo2 / (volt * denom + co2))
    return o2pct * o2sat_den / 100.0


def temp_rinko(temp_coeffs: Sequence[float], volt: float) -> float:
    # Temperature polynomial from the RINKO calibration coefficients.
    at, bt, ct, dt = temp_coeffs
    return np.poly1d([dt, ct, bt, at])(volt)


def elevation_corr(flag_pr: int, elevation_m: float, atmpr_mbar: float) -> float:
    # Air-pressure correction factor.
    if flag_pr == 0:
        return (1013.25 - elevation_m * 0.11568) / 1013.25
    return atmpr_mbar / 1013.25


def _time_factor(hour: float, time_points: Sequence[float], factors: Sequence[float]) -> float:
    # Piecewise-linear correction in time (lines 25-30 in the .def file).
    time1, time2, time3 = time_points
    f1, f2, f3 = factors
    if hour <= time2:
        denom = (time2 - time1)
        if denom == 0:
            return f2
        return f1 + (hour - time1) * (f2 - f1) / denom
    denom = (time3 - time2)
    if denom == 0:
        return f3
    return f2 + (hour - time2) * (f3 - f2) / denom


def _conc_factor(o2: float, o2_points: Sequence[float], factors: Sequence[float]) -> float:
    # Piecewise-linear correction based on oxygen concentration (lines 31-34).
    o24, o25 = o2_points
    f4, f5 = factors
    denom = (o25 - o24)
    if denom == 0:
        return f4
    return f4 + (o2 - o24) * (f5 - f4) / denom


def _iter_vector_rows(path: Path, column_map: ColumnMap) -> Iterable[Dict[str, float]]:
    # Stream rows from the converted Vector file to keep memory bounded.
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_num, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            parts = line.replace(",", " ").split()
            if parts and parts[0] == "\ufeff":
                parts = parts[1:]
            try:
                yield {
                    "vx": float(parts[column_map.vx - 1]),
                    "vy": float(parts[column_map.vy - 1]),
                    "vz": float(parts[column_map.vz - 1]),
                    "pr": float(parts[column_map.pr - 1]),
                    "p_counts": float(parts[column_map.p_counts - 1]),
                    "q_counts": float(parts[column_map.q_counts - 1]),
                }
            except IndexError as exc:
                raise ValueError(f"Line {line_num} has too few columns") from exc


def _resample_average(series: List[float], nav: int) -> List[float]:
    # Block-average to a lower frequency (Fortran-equivalent).
    out: List[float] = []
    n = (len(series) // nav) * nav
    for i in range(0, n, nav):
        out.append(sum(series[i:i + nav]) / nav)
    return out


def _resample_downsample(series: List[float], nav: int) -> List[float]:
    # Pick the first sample in each block.
    n = (len(series) // nav) * nav
    return [series[i] for i in range(0, n, nav)]


def _resample_interpolate(hours: List[float], series: List[float], target_hz: float) -> List[float]:
    # Linear interpolation onto a uniform target grid.
    if len(hours) < 2:
        return series[:]
    start = hours[0]
    end = hours[-1]
    total_seconds = (end - start) * 3600.0
    count = int(math.floor(total_seconds * target_hz)) + 1
    target_hours = [start + i / (target_hz * 3600.0) for i in range(count)]
    out: List[float] = []

    j = 0
    for t in target_hours:
        while j + 1 < len(hours) and hours[j + 1] < t:
            j += 1
        if t <= hours[0]:
            out.append(series[0])
        elif t >= hours[-1]:
            out.append(series[-1])
        else:
            t0 = hours[j]
            t1 = hours[j + 1]
            v0 = series[j]
            v1 = series[j + 1]
            if t1 == t0:
                out.append(v0)
            else:
                out.append(v0 + (v1 - v0) * (t - t0) / (t1 - t0))
    return out


def _resample_all(
    series: Dict[str, List[float]],
    hz: float,
    target_hz: float,
    mode: str,
) -> Dict[str, List[float]]:
    # Dispatch to the selected resampling strategy.
    if target_hz == hz:
        return series
    if target_hz > hz:
        raise ValueError("target_hz must be <= input hz")

    nav = int(round(hz / target_hz))
    if nav <= 0:
        raise ValueError("Invalid target_hz resulting in nav <= 0")

    out: Dict[str, List[float]] = {}
    if mode == "average":
        for name, values in series.items():
            out[name] = _resample_average(values, nav)
    elif mode == "downsample":
        for name, values in series.items():
            out[name] = _resample_downsample(values, nav)
    elif mode == "interpolate":
        hours = series["hour"]
        out["hour"] = _resample_interpolate(hours, hours, target_hz)
        for name, values in series.items():
            if name == "hour":
                continue
            out[name] = _resample_interpolate(hours, values, target_hz)
    else:
        raise ValueError(f"Unknown resample mode: {mode}")

    return out


def _format_row(values: Sequence[float]) -> str:
    # Single-line numeric output, space-separated.
    return " ".join(DEFAULT_FLOAT_FORMAT.format(v) for v in values)


def _chunk_output_path(base: Path, chunk_index: int) -> Path:
    # Output file per 30-minute (or configured) chunk.
    return base.with_name(f"{base.stem}_chunk{chunk_index:04d}{base.suffix}")


def process_definition(
    ini_path: Path,
    *,
    output_columns: Sequence[str] = OUTPUT_COLUMNS,
    chunk_minutes: int = 30,
    resample_mode: str = "average",
    target_hz: float | None = None,
    column_map: ColumnMap = ColumnMap(),
) -> None:
    # Main pipeline: read input, compute derived series, resample, and chunk outputs.
    definition = parse_ini_file(ini_path)

    input_path = Path(definition.input_file)
    if not input_path.is_absolute():
        input_path = ini_path.parent / input_path

    output_path = Path(definition.output_file)
    if not output_path.is_absolute():
        output_path = ini_path.parent / output_path

    if target_hz is None:
        # Default: no resampling unless explicitly requested.
        target_hz = definition.hz

    chunk_samples = int(round(definition.hz * (chunk_minutes * 60)))
    if chunk_samples <= 0:
        raise ValueError("chunk_minutes must be > 0")

    factor = elevation_corr(definition.flag_pr, definition.elevation_m, definition.atmpr_mbar)

    buffer: Dict[str, List[float]] = {
        "vx": [],
        "vy": [],
        "vz": [],
        "pr": [],
        "p_counts": [],
        "q_counts": [],
    }
    sample_index = 0
    chunk_index = 0

    def _series_value(series_map: Dict[str, List[float]], name: str, index: int) -> float:
        # Late validation for output columns.
        try:
            return series_map[name][index]
        except KeyError as exc:
            raise KeyError(f"Unknown output column: {name}") from exc

    def flush_chunk() -> None:
        nonlocal chunk_index
        if not buffer["vx"]:
            return

        hours: List[float] = []
        vx: List[float] = []
        vy: List[float] = []
        vz: List[float] = []
        pr_mbar: List[float] = []
        temp_c: List[float] = []
        o2sat_umol_l: List[float] = []
        o2_umol_l: List[float] = []

        start_index = sample_index - len(buffer["vx"])

        for i in range(len(buffer["vx"])):
            # Hour computed from start time and sample index.
            hour = definition.hour_start + ((start_index + i + 1) / definition.hz) / 3600.0
            hours.append(hour)

            vx.append(buffer["vx"][i] * 100.0)
            vy.append(buffer["vy"][i] * 100.0)
            vz.append(buffer["vz"][i] * 100.0)
            pr_mbar.append(buffer["pr"][i] * 100.0 - definition.pr_offset)

            volt1 = 5.0 * buffer["p_counts"][i] / 65535.0
            volt2 = 5.0 * buffer["q_counts"][i] / 65535.0

            temp = temp_rinko(definition.temp_coeffs, volt2)
            temp_c.append(temp)
            sat_den = o2sat(temp, definition.salinity) * density(temp, definition.salinity) * factor
            o2sat_umol_l.append(sat_den)

            o2_val = o2_rinko(
                definition.o2_coeffs,
                temp,
                volt1,
                sat_den,
            )
            o2_val *= _time_factor(
                hour,
                definition.time_points,
                definition.time_factors,
            )
            o2_val *= _conc_factor(o2_val, definition.o2_points, definition.o2_factors)
            o2_umol_l.append(o2_val)

        series: Dict[str, List[float]] = {
            "hour": hours,
            "vx": vx,
            "vy": vy,
            "vz": vz,
            "pr_mbar": pr_mbar,
            "temp_c": temp_c,
            "o2_umol_l": o2_umol_l,
            "o2sat_umol_l": o2sat_umol_l,
        }

        series = _resample_all(series, definition.hz, target_hz, resample_mode)

        out_path = _chunk_output_path(output_path, chunk_index)
        with out_path.open("w", encoding="utf-8", newline="\n") as handle:
            handle.write(" ".join(output_columns) + "\n")
            length = len(series["hour"])
            for i in range(length):
                row = [_series_value(series, col, i) for col in output_columns]
                handle.write(_format_row(row) + "\n")

        chunk_index += 1

    for row in _iter_vector_rows(input_path, column_map):
        buffer["vx"].append(row["vx"])
        buffer["vy"].append(row["vy"])
        buffer["vz"].append(row["vz"])
        buffer["pr"].append(row["pr"])
        buffer["p_counts"].append(row["p_counts"])
        buffer["q_counts"].append(row["q_counts"])
        sample_index += 1

        if len(buffer["vx"]) >= chunk_samples:
            flush_chunk()
            for key in buffer:
                buffer[key].clear()

    flush_chunk()


def main() -> int:
    parser = argparse.ArgumentParser(description="ExtractRINKO Python port")
    parser.add_argument("definition", type=Path, help="Path to the vector_formatter.ini file")
    parser.add_argument("--chunk-minutes", type=int, default=30, help="Chunk length in minutes")
    parser.add_argument(
        "--resample",
        choices=["downsample", "average", "interpolate"],
        default="average",
        help="Resampling mode for downscaling",
    )
    parser.add_argument("--target-hz", type=float, default=None, help="Target frequency (Hz).")
    resample_group = parser.add_mutually_exclusive_group()
    resample_group.add_argument(
        "--enable-resample",
        action="store_true",
        dest="enable_resample",
        help="Enable resampling.",
    )
    resample_group.add_argument(
        "--no-resample",
        action="store_false",
        dest="enable_resample",
        help="Do not resample (default).",
    )
    parser.set_defaults(enable_resample=False)
    args = parser.parse_args()

    if args.enable_resample:
        target_hz = args.target_hz
    else:
        target_hz = None

    process_definition(
        args.definition,
        chunk_minutes=args.chunk_minutes,
        resample_mode=args.resample,
        target_hz=target_hz,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
