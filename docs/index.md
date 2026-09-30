# wsu-sddata-generator

!!! note "Credits & Upstream Origins"
    This project is a fork of [wsu-sddata-generator](https://github.com/tnakazato/wsu-sddata-generator), created and authored by **Takeshi Nakazato ([@tnakazato](https://github.com/tnakazato))** at the National Astronomical Observatory of Japan (NAOJ). Full credit to **Takeshi Nakazato** for designing and implementing the original ALMA-WSU Single Dish data generator, including the MeasurementSet table manipulation architecture, channel expansion algorithms, and spectral window duplication logic.

    This fork extends the tool with vectorized visibility processing, modern Python 3.12+ packaging using `uv`, and automated documentation workflows.

`wsusd` generates synthetic ALMA Single Dish (SD) MeasurementSets emulating the Wideband Sensitivity Upgrade (WSU). Given an existing ASDM or MS, it modifies spectral resolution and bandwidth coverage.

## Overview

- **Channel expansion (`--chan-factor`)**: Increases the channel count for science spectral windows. Interpolates visibility data, scales weight and sigma spectra, interpolates `SYSCAL` calibration spectra, and injects residual noise.
- **SPW expansion (`--spw-factor`)**: Duplicates science and atmospheric spectral windows across MS sub-tables to simulate wider frequency coverage.
- **ASDM ingestion**: Automatically converts raw ASDM input using CASA `importasdm` before applying modifications.

## Documentation

- [Quickstart](quickstart.md): Environment setup and basic examples
- [Design](design.md): MS table manipulation details, interpolation, and noise modeling
- [CLI Reference](cli.md): Command-line options and syntax
