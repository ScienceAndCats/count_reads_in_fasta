# Count FASTQ reads against FASTA references with minimap2

This repository contains a Python command-line helper that aligns every FASTQ file in one folder against every FASTA file in another folder using [`minimap2`](https://github.com/lh3/minimap2), then writes a CSV summary of how many reads map to each reference after configurable quality filters.

## What the script does

`minimap2_fastq_vs_fasta_tally.py` performs an all-by-all comparison between input reads and reference sequences:

1. Reads settings from an INI configuration file.
2. Finds FASTQ files in `fastq_dir` with these extensions: `.fastq`, `.fq`, `.fastq.gz`, or `.fq.gz`.
3. Finds FASTA files in `fasta_dir` with these extensions: `.fasta`, `.fa`, `.fna`, `.fas`, and their `.gz` equivalents.
4. Counts total reads in each FASTQ and samples read lengths to choose a minimap2 preset when `preset = auto`.
5. Builds or reuses a minimap2 `.mmi` index for each FASTA/preset combination.
6. Aligns each FASTQ to each FASTA with minimap2 in SAM mode.
7. Counts primary alignments that pass the configured filters:
   - secondary and supplementary alignments are ignored;
   - unmapped reads are ignored;
   - alignments below `min_mapq` are ignored;
   - alignments without an `NM:i:` edit-distance tag are ignored;
   - alignments must have `NM / aligned_query_bases <= max_edit_fraction`;
   - alignments must have `aligned_query_bases / read_length >= min_query_covered`.
8. Writes one CSV row per FASTQ/FASTA pair.

The output CSV includes:

- `fastq_file`
- `fasta_file`
- `total_reads`
- `aligned_reads`
- `percent_aligned`
- `sampled_median_read_length`
- `preset_used`
- `max_edit_fraction`
- `min_query_covered`
- `min_mapq`

## Requirements

- Python 3.9 or newer is recommended.
- `minimap2` must be installed and available on your `PATH`.

The Python script only uses the Python standard library. No Python package installation is required.

### Install minimap2

Using conda/mamba:

```bash
mamba install -c bioconda minimap2
```

or:

```bash
conda install -c bioconda minimap2
```

Confirm that minimap2 is available:

```bash
minimap2 --version
```

## Repository files

- `minimap2_fastq_vs_fasta_tally.py` — the main script.
- `minimap2_tally_config.ini` — example/default configuration file used by the script.

## Quick start

1. Clone or download this repository.
2. Install `minimap2` and make sure it is on your `PATH`.
3. Edit `minimap2_tally_config.ini` so the paths point to your data.
4. Run the script:

```bash
python3 minimap2_fastq_vs_fasta_tally.py
```

By default, the script looks for `minimap2_tally_config.ini` in the same directory as `minimap2_fastq_vs_fasta_tally.py`.

## Configuration

Edit `minimap2_tally_config.ini` before running. The file has two sections: `[paths]` and `[settings]`.

### `[paths]`

```ini
[paths]
fastq_dir = /path/to/fastq_folder
fasta_dir = /path/to/fasta_folder
out_csv = minimap2_tally.csv
index_dir = minimap2_indexes
```

- `fastq_dir`: folder containing input FASTQ reads.
- `fasta_dir`: folder containing FASTA reference files.
- `out_csv`: path to the CSV summary that will be created.
- `index_dir`: folder where minimap2 `.mmi` indexes will be written and reused.

Relative paths are resolved relative to the configuration file location.

### `[settings]`

```ini
[settings]
preset = auto
long_read_preset = map-ont
recursive = false
threads = 4
max_edit_fraction = 0.10
min_query_covered = 0.00
min_mapq = 0
sample_reads_for_length = 1000
extra_mm2_args =
```

- `preset`: minimap2 preset to use. Allowed values are `auto`, `sr`, `map-ont`, `map-pb`, `asm5`, `asm10`, and `asm20`.
  - `auto` chooses `sr` when the sampled median read length is `<= 500` bases.
  - `auto` chooses `long_read_preset` when the sampled median read length is `> 500` bases.
- `long_read_preset`: long-read preset used only when `preset = auto`; choose `map-ont` or `map-pb`.
- `recursive`: set to `true` to search subfolders under `fastq_dir` and `fasta_dir`; otherwise only the top-level folders are searched.
- `threads`: number of CPU threads passed to minimap2.
- `max_edit_fraction`: maximum allowed edit fraction, calculated from minimap2's `NM` tag divided by aligned query bases. For example, `0.10` allows up to 10% mismatches/indels in the aligned part of the read.
- `min_query_covered`: minimum fraction of each read that must be aligned. For example, `0.80` requires at least 80% of the read to align.
- `min_mapq`: minimum MAPQ score required to count an alignment.
- `sample_reads_for_length`: number of reads sampled from each FASTQ to estimate median read length for automatic preset selection.
- `extra_mm2_args`: optional extra minimap2 arguments, exactly as you would type them on the command line, such as `-k 15 -w 10`.

## Running with a custom config file

You can keep multiple config files for different projects and choose one with `--config`:

```bash
python3 minimap2_fastq_vs_fasta_tally.py --config path/to/my_project_config.ini
```

## Example workflow

Suppose your project looks like this:

```text
project/
├── reads/
│   ├── sample_a.fastq.gz
│   └── sample_b.fastq.gz
├── references/
│   ├── target_1.fasta
│   └── target_2.fasta
└── config.ini
```

A matching `config.ini` could be:

```ini
[paths]
fastq_dir = reads
fasta_dir = references
out_csv = results/minimap2_tally.csv
index_dir = results/minimap2_indexes

[settings]
preset = auto
long_read_preset = map-ont
recursive = false
threads = 8
max_edit_fraction = 0.10
min_query_covered = 0.80
min_mapq = 10
sample_reads_for_length = 1000
extra_mm2_args =
```

Run from the `project/` directory:

```bash
python3 /path/to/count_reads_in_fasta/minimap2_fastq_vs_fasta_tally.py --config config.ini
```

The script will create `results/minimap2_tally.csv` with one row for each FASTQ/FASTA combination.

## Notes and limitations

- The script counts primary SAM records only. Secondary (`0x100`) and supplementary (`0x800`) alignments are skipped.
- A read is counted for a FASTA only if its primary alignment passes all filters.
- The script compares every FASTQ against every FASTA. Runtime scales with `number of FASTQ files × number of FASTA files`.
- Existing `.mmi` indexes are reused unless the FASTA file has been modified more recently than the index.
- If you change minimap2 arguments that affect indexing, delete old indexes from `index_dir` before rerunning.
- Input FASTQ and FASTA files may be plain text or gzip-compressed.

## Troubleshooting

### `minimap2 was not found in PATH`

Install minimap2 and confirm that `minimap2 --version` works in the same terminal where you run the Python script.

### `FASTQ directory not found` or `FASTA directory not found`

Check the paths in `[paths]`. Relative paths are resolved relative to the config file, not necessarily relative to your current shell directory.

### No FASTQ or FASTA files found

Confirm that your files use one of the supported extensions and that `recursive = true` if your files are inside subfolders.

### Results are too strict or too permissive

Adjust `max_edit_fraction`, `min_query_covered`, and `min_mapq`. For stricter matching, lower `max_edit_fraction` and raise `min_query_covered` or `min_mapq`.
