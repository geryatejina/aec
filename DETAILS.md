# ExtractRINKO Processing Notes

This file captures the essential input/output structure and processing flow we need to preserve while porting ExtractRINKO from Fortran to Python. It is intentionally more detailed than the README.

## Inputs
ExtractRINKO uses two inputs:
1) **Converted Nortek Vector data file** (text, 18 columns).
2) **Definition file** (text) that controls processing and output formatting.

### Vector data file (converted from binary)
Key columns (1-based, as described in the manual):
- 3–5: velocity components (x, y, z)
- 15: pressure
- 16: oxygen signal (counts)
- 17: temperature signal (counts)

### Definition file (line-based controls)
Important lines (1-based):
- 2: Vector data filename to process
- 3: Output filename
- 4: Start time (hours)
- 5: Input sampling frequency (Hz)
- 6: Output sampling frequency (Hz)
- 7: Pressure offset (mbar)
- 8: Salinity (ppt)
- 9–11: Atmospheric pressure handling (either elevation or measured pressure)
- 12–15: Temperature conversion coefficients (sensor-specific)
- 16–23: Oxygen conversion coefficients (sensor/foil-specific)
- 24: Output format flag (controls columns 5–6)
- 25–30: Optional time-varying calibration factor (piecewise linear)
- 31–34: Optional concentration-based calibration factor (linear)

## Output file
Output always contains 6 columns:
1) time (hours)
2–4) velocity (x, y, z) in cm/s
5–6) depend on the output flag in line 24 of the definition file (raw counts or computed values)

## Processing overview (high-level)
- Convert counts to voltage.
- Compute high-frequency temperature from coefficients.
- Compute oxygen saturation and oxygen concentration, corrected by temperature, salinity, and air pressure.
- Apply optional calibration factors (time-based or concentration-based).

## Porting checklist (draft)
- Parse definition file with strict validation and clear errors.
- Read Vector file with explicit column handling and unit conversions.
- Reproduce numerical formulas; add unit tests and regression fixtures.
- Provide a CLI and a Python API.

## Notes
- The Python port should preserve numerical results within a tight tolerance.
- All unit conversions and constants must be documented and tested.
