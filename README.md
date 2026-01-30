# ExtractRINKO (Python Port)

This repo will port the ExtractRINKO Fortran workflow to Python and add a few quality-of-life improvements (validation, clearer I/O, tests, and reproducible processing).

## Goals
- Reproduce ExtractRINKO results from the Nortek Vector + RINKO EC sensor pipeline.
- Provide a clean, well-tested Python implementation.
- Make inputs/outputs explicit and easier to validate.

## Status
Early setup: repo scaffolding and documentation.

## Project layout (planned)
- `src/`: Python implementation
- `docs/`: reference documents (ignored by git)
- `sample_data/`: example inputs/outputs (external)

## Quick start (placeholder)
Configuration now uses `vector_formatter.ini`. You can convert legacy `.def` files with:
```
python src/def_to_ini.py path\to\file.def
```

Basic usage (default = no resampling):
```
python src/vector_formatter.py path\to\vector_formatter.ini
```

Enable resampling (downsample/average/interpolate):
```
python src/vector_formatter.py path\to\vector_formatter.ini --enable-resample --resample average --target-hz 16
```

## References
See `DETAILS.md` for a concise breakdown of inputs/outputs and processing steps.
