#!/usr/bin/env python3
from __future__ import annotations

import argparse
import configparser
import csv
import gzip
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from statistics import median


FASTQ_EXTENSIONS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
FASTA_EXTENSIONS = (".fasta", ".fa", ".fna", ".fas", ".fasta.gz", ".fa.gz", ".fna.gz", ".fas.gz")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
DEFAULT_CONFIG_NAME = "minimap2_tally_config.ini"


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def find_files(folder: Path, extensions: tuple[str, ...], recursive: bool = False) -> list[Path]:
    files: list[Path] = []
    iterator = folder.rglob("*") if recursive else folder.iterdir()
    for p in iterator:
        if p.is_file() and any(str(p).lower().endswith(ext) for ext in extensions):
            files.append(p)
    return sorted(files)


def fastq_stats(path: Path, sample_reads: int = 1000) -> tuple[int, float]:
    total_reads = 0
    lengths: list[int] = []
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                raise ValueError(f"Malformed FASTQ file: {path}")
            total_reads += 1
            if len(lengths) < sample_reads:
                lengths.append(len(seq.rstrip("\n\r")))
    med = float(median(lengths)) if lengths else 0.0
    return total_reads, med


def choose_preset(user_preset: str, median_len: float, long_read_preset: str) -> str:
    if user_preset != "auto":
        return user_preset
    return "sr" if median_len <= 500 else long_read_preset


def ensure_minimap2() -> str:
    exe = shutil.which("minimap2")
    if exe is None:
        raise FileNotFoundError(
            "minimap2 was not found in PATH. Install it first, e.g. conda install -c bioconda minimap2"
        )
    return exe


def nm_from_fields(fields: list[str]) -> int | None:
    for field in fields[11:]:
        if field.startswith("NM:i:"):
            return int(field[5:])
    return None


def aligned_query_bases(cigar: str) -> int:
    total = 0
    for n_str, op in CIGAR_RE.findall(cigar):
        n = int(n_str)
        if op in {"M", "I", "=", "X"}:
            total += n
    return total


def read_length_from_cigar_and_seq(cigar: str, seq: str) -> int:
    if seq != "*":
        return len(seq)
    total = 0
    for n_str, op in CIGAR_RE.findall(cigar):
        n = int(n_str)
        if op in {"M", "I", "S", "=", "X"}:
            total += n
    return total


def ensure_index(minimap2_exe: str, fasta_path: Path, index_dir: Path, preset: str) -> Path:
    index_dir.mkdir(parents=True, exist_ok=True)
    index_path = index_dir / f"{fasta_path.name}.{preset}.mmi"
    needs_build = (not index_path.exists()) or (index_path.stat().st_mtime < fasta_path.stat().st_mtime)
    if needs_build:
        cmd = [minimap2_exe, "-x", preset, "-d", str(index_path), str(fasta_path)]
        subprocess.run(cmd, check=True)
    return index_path


def align_and_tally(
    minimap2_exe: str,
    fastq_path: Path,
    fasta_path: Path,
    index_path: Path,
    preset: str,
    threads: int,
    max_edit_fraction: float,
    min_query_covered: float,
    min_mapq: int,
    extra_args: list[str],
) -> tuple[int, int]:
    cmd = [
        minimap2_exe,
        "-a",
        "-x",
        preset,
        "-t",
        str(threads),
        *extra_args,
        str(index_path),
        str(fastq_path),
    ]

    total_primary = 0
    passing_primary = 0

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    assert proc.stdout is not None

    for line in proc.stdout:
        if not line or line.startswith("@"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 11:
            continue

        flag = int(fields[1])
        if flag & 0x100 or flag & 0x800:
            continue

        total_primary += 1

        if flag & 0x4:
            continue

        mapq = int(fields[4])
        if mapq < min_mapq:
            continue

        cigar = fields[5]
        seq = fields[9]
        nm = nm_from_fields(fields)
        if nm is None:
            continue

        q_aln_bases = aligned_query_bases(cigar)
        q_len = read_length_from_cigar_and_seq(cigar, seq)
        if q_aln_bases <= 0 or q_len <= 0:
            continue

        edit_fraction = nm / q_aln_bases
        query_covered = q_aln_bases / q_len

        if edit_fraction <= max_edit_fraction and query_covered >= min_query_covered:
            passing_primary += 1

    stderr_text = proc.stderr.read() if proc.stderr is not None else ""
    return_code = proc.wait()
    if return_code != 0:
        raise RuntimeError(
            f"minimap2 failed for FASTQ={fastq_path.name} vs FASTA={fasta_path.name}\n"
            f"Command: {' '.join(cmd)}\n\n{stderr_text}"
        )

    return total_primary, passing_primary


def resolve_path(config_dir: Path, value: str) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = (config_dir / path).resolve()
    return path


def load_config(config_path: Path) -> dict[str, object]:
    parser = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    read_files = parser.read(config_path)
    if not read_files:
        raise FileNotFoundError(f"Could not read config file: {config_path}")

    required_sections = {"paths", "settings"}
    missing_sections = [s for s in required_sections if s not in parser]
    if missing_sections:
        raise ValueError(f"Config file is missing required section(s): {', '.join(missing_sections)}")

    cfg_dir = config_path.resolve().parent

    def get_required(section: str, option: str) -> str:
        if option not in parser[section] or not parser[section][option].strip():
            raise ValueError(f"Config file is missing required value: [{section}] {option}")
        return parser[section][option].strip()

    fastq_dir = resolve_path(cfg_dir, get_required("paths", "fastq_dir"))
    fasta_dir = resolve_path(cfg_dir, get_required("paths", "fasta_dir"))
    out_csv = resolve_path(cfg_dir, get_required("paths", "out_csv"))
    index_dir = resolve_path(cfg_dir, get_required("paths", "index_dir"))

    preset = parser.get("settings", "preset", fallback="auto").strip()
    long_read_preset = parser.get("settings", "long_read_preset", fallback="map-ont").strip()
    recursive = parser.getboolean("settings", "recursive", fallback=False)
    threads = parser.getint("settings", "threads", fallback=4)
    max_edit_fraction = parser.getfloat("settings", "max_edit_fraction", fallback=0.10)
    min_query_covered = parser.getfloat("settings", "min_query_covered", fallback=0.0)
    min_mapq = parser.getint("settings", "min_mapq", fallback=0)
    sample_reads_for_length = parser.getint("settings", "sample_reads_for_length", fallback=1000)
    extra_mm2_args_raw = parser.get("settings", "extra_mm2_args", fallback="").strip()
    extra_mm2_args = shlex.split(extra_mm2_args_raw) if extra_mm2_args_raw else []

    allowed_presets = {"auto", "sr", "map-ont", "map-pb", "asm5", "asm10", "asm20"}
    if preset not in allowed_presets:
        raise ValueError(f"Invalid preset '{preset}'. Allowed: {', '.join(sorted(allowed_presets))}")
    if long_read_preset not in {"map-ont", "map-pb"}:
        raise ValueError("long_read_preset must be 'map-ont' or 'map-pb'")
    if not (0.0 <= max_edit_fraction <= 1.0):
        raise ValueError("max_edit_fraction must be between 0 and 1")
    if not (0.0 <= min_query_covered <= 1.0):
        raise ValueError("min_query_covered must be between 0 and 1")
    if min_mapq < 0:
        raise ValueError("min_mapq must be >= 0")
    if threads < 1:
        raise ValueError("threads must be >= 1")
    if sample_reads_for_length < 1:
        raise ValueError("sample_reads_for_length must be >= 1")

    return {
        "config_path": config_path.resolve(),
        "fastq_dir": fastq_dir,
        "fasta_dir": fasta_dir,
        "out_csv": out_csv,
        "index_dir": index_dir,
        "preset": preset,
        "long_read_preset": long_read_preset,
        "recursive": recursive,
        "threads": threads,
        "max_edit_fraction": max_edit_fraction,
        "min_query_covered": min_query_covered,
        "min_mapq": min_mapq,
        "sample_reads_for_length": sample_reads_for_length,
        "extra_mm2_args": extra_mm2_args,
    }


def get_config_path() -> Path:
    script_dir = Path(__file__).resolve().parent
    default_config = script_dir / DEFAULT_CONFIG_NAME

    ap = argparse.ArgumentParser(
        description=(
            "Align every FASTQ in a folder against every FASTA in another folder with minimap2, "
            "then count how many primary reads align to each FASTA after filtering by edit fraction. "
            "Edit the config file instead of typing long command lines."
        )
    )
    ap.add_argument(
        "--config",
        default=str(default_config),
        help=(
            "Path to config INI file. Default: a file named "
            f"'{DEFAULT_CONFIG_NAME}' next to this script."
        ),
    )
    args = ap.parse_args()
    return Path(args.config).expanduser().resolve()


def main() -> None:
    config_path = get_config_path()
    config = load_config(config_path)

    minimap2_exe = ensure_minimap2()
    fastq_dir = config["fastq_dir"]
    fasta_dir = config["fasta_dir"]
    out_csv = config["out_csv"]
    index_dir = config["index_dir"]
    preset_setting = config["preset"]
    long_read_preset = config["long_read_preset"]
    recursive = config["recursive"]
    threads = config["threads"]
    max_edit_fraction = config["max_edit_fraction"]
    min_query_covered = config["min_query_covered"]
    min_mapq = config["min_mapq"]
    sample_reads_for_length = config["sample_reads_for_length"]
    extra_mm2_args = config["extra_mm2_args"]

    assert isinstance(fastq_dir, Path)
    assert isinstance(fasta_dir, Path)
    assert isinstance(out_csv, Path)
    assert isinstance(index_dir, Path)

    if not fastq_dir.exists():
        raise FileNotFoundError(f"FASTQ directory not found: {fastq_dir}")
    if not fasta_dir.exists():
        raise FileNotFoundError(f"FASTA directory not found: {fasta_dir}")

    fastq_files = find_files(fastq_dir, FASTQ_EXTENSIONS, recursive=bool(recursive))
    fasta_files = find_files(fasta_dir, FASTA_EXTENSIONS, recursive=bool(recursive))

    if not fastq_files:
        raise FileNotFoundError(f"No FASTQ files found in {fastq_dir}")
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in {fasta_dir}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []

    print(f"Config: {config_path}", flush=True)
    print(f"Found {len(fastq_files)} FASTQ file(s) and {len(fasta_files)} FASTA file(s).", flush=True)
    print(
        "Counting primary mapped reads only; secondary and supplementary alignments are ignored.",
        flush=True,
    )
    print(
        f"Filters: max_edit_fraction={max_edit_fraction}, min_query_covered={min_query_covered}, min_mapq={min_mapq}",
        flush=True,
    )
    print(flush=True)

    for fastq_path in fastq_files:
        total_reads, median_len = fastq_stats(fastq_path, sample_reads=int(sample_reads_for_length))
        preset = choose_preset(str(preset_setting), float(median_len), str(long_read_preset))

        print(
            f"FASTQ: {fastq_path.name} | total_reads={total_reads} | "
            f"sampled_median_read_len={median_len:.1f} | preset={preset}",
            flush=True,
        )

        for fasta_path in fasta_files:
            index_path = ensure_index(minimap2_exe, fasta_path, index_dir, preset)
            total_primary, aligned_reads = align_and_tally(
                minimap2_exe=minimap2_exe,
                fastq_path=fastq_path,
                fasta_path=fasta_path,
                index_path=index_path,
                preset=preset,
                threads=int(threads),
                max_edit_fraction=float(max_edit_fraction),
                min_query_covered=float(min_query_covered),
                min_mapq=int(min_mapq),
                extra_args=list(extra_mm2_args),
            )

            if total_primary != total_reads:
                print(
                    f"  WARNING: primary SAM records ({total_primary}) != FASTQ reads ({total_reads}) for {fastq_path.name}",
                    file=sys.stderr,
                    flush=True,
                )

            percent_aligned = (aligned_reads / total_reads * 100.0) if total_reads else 0.0
            print(
                f"  vs {fasta_path.name}: aligned_reads={aligned_reads}/{total_reads} ({percent_aligned:.2f}%)",
                flush=True,
            )

            rows.append(
                {
                    "fastq_file": fastq_path.name,
                    "fasta_file": fasta_path.name,
                    "total_reads": total_reads,
                    "aligned_reads": aligned_reads,
                    "percent_aligned": f"{percent_aligned:.4f}",
                    "sampled_median_read_length": f"{median_len:.1f}",
                    "preset_used": preset,
                    "max_edit_fraction": max_edit_fraction,
                    "min_query_covered": min_query_covered,
                    "min_mapq": min_mapq,
                }
            )
        print(flush=True)

    with open(out_csv, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "fastq_file",
                "fasta_file",
                "total_reads",
                "aligned_reads",
                "percent_aligned",
                "sampled_median_read_length",
                "preset_used",
                "max_edit_fraction",
                "min_query_covered",
                "min_mapq",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote CSV summary to: {out_csv}", flush=True)


if __name__ == "__main__":
    main()
