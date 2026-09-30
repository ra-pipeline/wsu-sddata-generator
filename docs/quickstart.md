# Quickstart

## Requirements

- Python >= 3.12
- [uv](https://docs.astral.sh/uv/) (recommended)
- CASA package index: `https://casa-pip.nrao.edu/repository/pypi-group/simple`

## Setup

The repository is configured for `uv`, with `casatools` and `casatasks` mapped directly to the CASA package index in `pyproject.toml`.

To create and populate the local `.venv`:

```bash
uv sync --all-groups
source .venv/bin/activate
```

Or run commands directly via `uv`:

```bash
uv run wsusd --help
```

## Examples

### Expand channels on an MS
Expand science spectral windows by a factor of 4:

```bash
wsusd --chan-factor 4 my_target.ms
```

To create a `.bak` copy of the MS before in-place modification:

```bash
wsusd -b -c 4 my_target.ms
```

### Process an ASDM directly
When passed an ASDM directory, `wsusd` runs CASA `importasdm` and names the output MS based on the expansion factors:

```bash
wsusd -c 2 -s 2 /path/to/uid___A002_X117d38c_X13867
# Writes: uid___A002_X117d38c_X13867.2.0xchan.2xspw.ms
```

### Dry run
Print parsed parameters and target SPW counts without modifying data:

```bash
wsusd --dry-run -c 4 -s 2 my_target.ms
```
