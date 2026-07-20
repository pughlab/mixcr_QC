#!/usr/bin/env python3
"""
Parse MiXCR LOG_ALIGN / LOG_ASSEMBLE files across multiple samples.

Expects the per-sample folder layout produced by run_batch_pipeline.sh /
submit_pipeline_slurm.sh:

    BASE_DIR/
        <sample_1>/LOG/LOG_ALIGN_<sample_1>.txt
        <sample_1>/LOG/LOG_ASSEMBLE_<sample_1>.txt
        <sample_2>/LOG/LOG_ALIGN_<sample_2>.txt
        <sample_2>/LOG/LOG_ASSEMBLE_<sample_2>.txt
        ...

Writes three files to --out_dir:
    align_stats.csv     - one row per sample, align-log columns (same
                           schema as the original log_parser.py output,
                           so qc_plotter.py can consume it directly)
    assemble_stats.csv  - one row per sample, assemble-log columns
    qc_metrics.tsv       - align + assemble merged into one row per sample,
                           for quick inspection

Parsing logic (key:value lines, chain zero-fill, Version/Input/Output file
handling, percentage stripping) is carried over from the original
single-directory log_parser.py, adapted to run per sample and merged rather
than relying on unjoined threads.

USAGE:
    python3 parse_and_merge_qc.py --base_dir OUT_DIR --out_dir OUT_DIR
"""

import argparse
import csv
import glob
import os


ALIGN_SKIP_PREFIXES = (
    "Analysis Date", "Overlapped", "Analysis time", "Command line arguments",
    "Paired-end alignment conflicts eliminated", "J gene chimeras",
    "V gene chimeras", "Chimeras",
)
ASSEMBLE_SKIP_PREFIXES = ("Analysis Date", "Analysis time", "Command line arguments")

ALIGN_CHAINS = ['TRA chains', 'TRB chains', 'TRD chains', 'TRG chains', 'TRA,TRD chains',
                'IGH chains', 'IGK chains', 'IGL chains']
ASSEMBLE_CHAINS = ["IGH chains", "IGK chains", "IGL chains",
                    "TRA chains", "TRB chains", "TRD chains", "TRG chains"]


def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base_dir", type=str, required=True,
                        help="base pipeline output directory containing one subfolder per sample")
    parser.add_argument("--out_dir", type=str, required=True,
                        help="directory to write align_stats.csv, assemble_stats.csv, qc_metrics.tsv")
    parser.add_argument("-pa", "--prefix_align", type=str, default="LOG_ALIGN",
                        help="prefix of alignment log file")
    parser.add_argument("-pm", "--prefix_assemble", type=str, default="LOG_ASSEMBLE",
                        help="prefix of assemble log file")
    return parser.parse_args()


def _parse_log_block(log_file, skip_prefixes, chain_list):
    """Parse a single LOG_ALIGN or LOG_ASSEMBLE file into one flat dict."""
    row = {}
    with open(log_file) as f:
        for line in f:
            try:
                if line.startswith(skip_prefixes):
                    continue
                if line.startswith("==="):
                    for chain in chain_list:
                        row.setdefault(chain, '0')
                    row.setdefault("chains", '0')
                    continue
                if "Version" not in line:
                    parts = line.strip().split(":")
                    if len(parts) < 2:
                        continue
                    key, val = parts[0], parts[1]
                else:
                    parts = line.strip().split(";")
                    key, val = "Version", parts[0].split(":")[1] + "," + parts[-1].split("=")[1]

                if key == 'Input file(s)':
                    val = ",".join(os.path.basename(v) for v in val.split(","))
                elif key == 'Output file':
                    val = os.path.basename(val)
                elif key == 'Version':
                    val = val.split(";")[0]
                elif "%)" in val:
                    val = val.split()[0]
                row[key] = val
            except Exception as e:
                print(f"error parsing line in {log_file}: {e}")
                continue
    return row


def parse_sample_logs(sample_log_dir, sample_name, prefix_align, prefix_assemble):
    """Parse the align + assemble logs for one sample; return (align_row, assemble_row)."""
    align_files = glob.glob(os.path.join(sample_log_dir, f"{prefix_align}*.txt"))
    assemble_files = glob.glob(os.path.join(sample_log_dir, f"{prefix_assemble}*.txt"))

    align_row = None
    assemble_row = None

    if align_files:
        align_row = _parse_log_block(align_files[0], ALIGN_SKIP_PREFIXES, ALIGN_CHAINS)
        align_row["Sample"] = sample_name
    else:
        print(f"WARNING: no {prefix_align}*.txt found for {sample_name}")

    if assemble_files:
        assemble_row = _parse_log_block(assemble_files[0], ASSEMBLE_SKIP_PREFIXES, ASSEMBLE_CHAINS)
        assemble_row["Sample"] = sample_name
    else:
        print(f"WARNING: no {prefix_assemble}*.txt found for {sample_name}")

    return align_row, assemble_row


def write_csv(rows, output_file, delimiter=','):
    """Write rows (list of dicts) to a delimited file, unioning columns across rows."""
    columns = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)

    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter=delimiter, restval='NA')
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def write_merged_tsv(align_rows, assemble_rows, output_file):
    """One row per sample, align_/assemble_-prefixed columns, for quick inspection."""
    by_sample = {}
    for row in align_rows:
        s = row["Sample"]
        by_sample.setdefault(s, {"Sample": s})
        by_sample[s].update({f"align_{k}": v for k, v in row.items() if k != "Sample"})
    for row in assemble_rows:
        s = row["Sample"]
        by_sample.setdefault(s, {"Sample": s})
        by_sample[s].update({f"assemble_{k}": v for k, v in row.items() if k != "Sample"})

    write_csv(list(by_sample.values()), output_file, delimiter='\t')


def main():
    args = get_options()
    os.makedirs(args.out_dir, exist_ok=True)

    sample_dirs = sorted(
        d for d in glob.glob(os.path.join(args.base_dir, "*"))
        if os.path.isdir(d) and os.path.isdir(os.path.join(d, "LOG"))
    )

    if not sample_dirs:
        raise SystemExit(f"No sample subfolders with a LOG/ directory found under {args.base_dir}")

    align_rows, assemble_rows = [], []
    for sample_dir in sample_dirs:
        sample_name = os.path.basename(sample_dir.rstrip("/"))
        log_dir = os.path.join(sample_dir, "LOG")
        print(f"Parsing logs for {sample_name}...")
        align_row, assemble_row = parse_sample_logs(log_dir, sample_name, args.prefix_align, args.prefix_assemble)
        if align_row:
            align_rows.append(align_row)
        if assemble_row:
            assemble_rows.append(assemble_row)

    align_stats_path = os.path.join(args.out_dir, "align_stats.csv")
    assemble_stats_path = os.path.join(args.out_dir, "assemble_stats.csv")
    merged_path = os.path.join(args.out_dir, "qc_metrics.tsv")

    write_csv(align_rows, align_stats_path)
    write_csv(assemble_rows, assemble_stats_path)
    write_merged_tsv(align_rows, assemble_rows, merged_path)

    print(f"Wrote {len(align_rows)} align row(s) to {align_stats_path}")
    print(f"Wrote {len(assemble_rows)} assemble row(s) to {assemble_stats_path}")
    print(f"Wrote merged QC metrics to {merged_path}")


if __name__ == "__main__":
    main()
