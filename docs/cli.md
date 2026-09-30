# CLI Reference

## Synopsis

```bash
wsusd [-h] [--version] [--debug] [--dry-run] [--backup-ms]
      [--chan-factor CHAN_FACTOR] [--spw-factor SPW_FACTOR]
      asdm_name
```

## Arguments

| Argument | Type | Description |
| :--- | :--- | :--- |
| `asdm_name` | Positional | Path to the input ASDM directory or MeasurementSet. |

## Options

| Option | Default | Description |
| :--- | :--- | :--- |
| `-c, --chan-factor` | `1` | Channel expansion factor (float > 0). Setting this to `10` produces 10x more channels than input. |
| `-s, --spw-factor` | `1` | Spectral window expansion factor (integer >= 1). Setting this to `2` duplicates science SPWs once. |
| `-b, --backup-ms` | `False` | Back up MS to `<vis>.bak` before modifying. |
| `--dry-run` | `False` | Print parsed input parameters and target output specs, then exit without modifying data. |
| `-d, --debug` | `False` | Enable debug logging. |
| `--version` | | Show version number and exit. |
| `-h, --help` | | Show help message and exit. |

## Examples

```bash
# Dry run on an ASDM
wsusd --dry-run -c 4 -s 2 /path/to/uid___A002_X117d38c_X13867

# In-place MS channel expansion with backup
wsusd -b -c 10 my_data.ms

# Run with debug logging
wsusd -d -c 2 my_data.ms
```
