#!/usr/bin/env python3

import sys
import csv
import argparse
from collections import Counter

DEFAULT_MISSING = {"N", "n", "-", "?"}
IUPAC_AMBIG = set("RYSWKMBDHVryswkmbdhv")


def read_fasta(path):
    sequences = []
    name = None
    seq_chunks = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    sequences.append((name, "".join(seq_chunks)))
                name = line[1:].strip()
                seq_chunks = []
            else:
                seq_chunks.append(line)

        if name is not None:
            sequences.append((name, "".join(seq_chunks)))

    return sequences


def write_fasta(path, header, sequence, line_width=80):
    with open(path, "w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i + line_width] + "\n")


def is_missing(char, treat_iupac_as_missing=False):
    if char in DEFAULT_MISSING:
        return True
    if treat_iupac_as_missing and char in IUPAC_AMBIG:
        return True
    return False


def build_priority_consensus(
    sequences,
    treat_iupac_as_missing=False,
    ignore_insertions_relative_to_first=True
):
    if not sequences:
        raise ValueError("No sequences found in FASTA file.")

    lengths = {len(seq) for _, seq in sequences}
    if len(lengths) != 1:
        raise ValueError("Sequences are not all the same length. Input must be an aligned FASTA.")

    aln_len = lengths.pop()
    nseq = len(sequences)

    first_name, first_seq = sequences[0]

    consensus = []
    contributors = []
    position_rows = []

    contribution_counts = Counter()
    data_counts_anchor_only = Counter()
    missing_counts_anchor_only = Counter()
    ambiguity_counts_anchor_only = Counter()
    base_counts_in_consensus = Counter()

    skipped_insertion_columns = 0
    total_anchor_positions = 0

    # Stats only across retained columns (i.e. columns kept in the final sequence)
    for pos in range(aln_len):
        if ignore_insertions_relative_to_first and first_seq[pos] == "-":
            skipped_insertion_columns += 1
            continue

        total_anchor_positions += 1

        for idx, (_, seq) in enumerate(sequences, start=1):
            c = seq[pos]
            if is_missing(c, treat_iupac_as_missing=treat_iupac_as_missing):
                missing_counts_anchor_only[idx] += 1
            else:
                data_counts_anchor_only[idx] += 1
                if c in IUPAC_AMBIG:
                    ambiguity_counts_anchor_only[idx] += 1

    anchor_pos_1based = 0

    for aln_pos in range(aln_len):
        # Ignore columns that are insertions relative to sequence 1
        if ignore_insertions_relative_to_first and first_seq[aln_pos] == "-":
            continue

        anchor_pos_1based += 1
        chosen_base = "N"
        chosen_seq_idx = 0
        bases_at_pos = [seq[aln_pos] for _, seq in sequences]

        for idx, (_, seq) in enumerate(sequences, start=1):
            c = seq[aln_pos]

            # never allow lower-priority inserted bases to create sequence
            # because columns with seq1='-' were already skipped above
            if not is_missing(c, treat_iupac_as_missing=treat_iupac_as_missing):
                chosen_base = c.upper()
                chosen_seq_idx = idx
                contribution_counts[idx] += 1
                break

        if chosen_seq_idx == 0:
            contribution_counts[0] += 1

        consensus.append(chosen_base)
        contributors.append(chosen_seq_idx)
        base_counts_in_consensus[chosen_base] += 1

        row = {
            "alignment_position_1based": aln_pos + 1,
            "output_position_1based": anchor_pos_1based,
            "consensus_base": chosen_base,
            "contributor_rank": chosen_seq_idx if chosen_seq_idx != 0 else "none",
            "contributor_name": sequences[chosen_seq_idx - 1][0] if chosen_seq_idx != 0 else "none",
            "seq1_anchor_base": first_seq[aln_pos],
        }

        for idx, base in enumerate(bases_at_pos, start=1):
            row[f"seq{idx}_base"] = base

        position_rows.append(row)

    return {
        "consensus": "".join(consensus),
        "contributors": contributors,
        "position_rows": position_rows,
        "contribution_counts": contribution_counts,
        "data_counts_anchor_only": data_counts_anchor_only,
        "missing_counts_anchor_only": missing_counts_anchor_only,
        "ambiguity_counts_anchor_only": ambiguity_counts_anchor_only,
        "base_counts_in_consensus": base_counts_in_consensus,
        "alignment_length": aln_len,
        "output_length": len(consensus),
        "n_sequences": nseq,
        "skipped_insertion_columns": skipped_insertion_columns,
        "anchor_name": first_name,
    }


def write_position_report(path, position_rows, n_sequences):
    fieldnames = [
        "alignment_position_1based",
        "output_position_1based",
        "consensus_base",
        "contributor_rank",
        "contributor_name",
        "seq1_anchor_base",
    ] + [f"seq{i}_base" for i in range(1, n_sequences + 1)]

    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in position_rows:
            writer.writerow(row)


def write_summary_report(path, sequences, results, treat_iupac_as_missing=False, ignore_insertions_relative_to_first=True):
    output_len = results["output_length"]
    aln_len = results["alignment_length"]
    n_sequences = results["n_sequences"]
    contribution_counts = results["contribution_counts"]
    data_counts = results["data_counts_anchor_only"]
    missing_counts = results["missing_counts_anchor_only"]
    ambiguity_counts = results["ambiguity_counts_anchor_only"]
    base_counts = results["base_counts_in_consensus"]
    skipped_insertion_columns = results["skipped_insertion_columns"]
    anchor_name = results["anchor_name"]

    total_filled = output_len - contribution_counts.get(0, 0)
    total_missing_in_consensus = contribution_counts.get(0, 0)

    with open(path, "w") as f:
        f.write("PRIORITY CONSENSUS SUMMARY REPORT\n")
        f.write("================================\n\n")

        f.write(f"Number of input sequences: {n_sequences}\n")
        f.write(f"Alignment length: {aln_len}\n")
        f.write(f"Anchor sequence (highest priority): {anchor_name}\n")
        f.write(f"Ignore insertions relative to first sequence: {'yes' if ignore_insertions_relative_to_first else 'no'}\n")
        f.write(f"Columns skipped because seq1 had '-': {skipped_insertion_columns}\n")
        f.write(f"Output consensus length: {output_len}\n")
        f.write(f"Consensus positions filled by at least one sequence: {total_filled}\n")
        f.write(f"Consensus positions remaining missing (all sequences missing): {total_missing_in_consensus}\n")
        f.write(f"Percent consensus resolved: {100 * total_filled / output_len:.2f}%\n" if output_len else "Percent consensus resolved: 0.00%\n")
        f.write(f"Percent consensus unresolved: {100 * total_missing_in_consensus / output_len:.2f}%\n" if output_len else "Percent consensus unresolved: 0.00%\n")
        f.write(f"Treat IUPAC ambiguity as missing: {'yes' if treat_iupac_as_missing else 'no'}\n\n")

        f.write("Consensus base composition\n")
        f.write("-------------------------\n")
        for base in ["A", "C", "G", "T", "N"]:
            count = base_counts.get(base, 0)
            pct = 100 * count / output_len if output_len else 0
            f.write(f"{base}: {count} ({pct:.2f}%)\n")

        other_bases = sorted(b for b in base_counts.keys() if b not in {"A", "C", "G", "T", "N"})
        for base in other_bases:
            count = base_counts[base]
            pct = 100 * count / output_len if output_len else 0
            f.write(f"{base}: {count} ({pct:.2f}%)\n")

        f.write("\nPer-sequence contribution\n")
        f.write("-------------------------\n")
        f.write(
            "rank\tname\tcontributed_positions\t%_of_output\t%_of_resolved_consensus\t"
            "positions_with_data_in_retained_columns\tpositions_missing_in_retained_columns\tambiguity_codes_in_retained_columns\n"
        )

        for idx, (name, seq) in enumerate(sequences, start=1):
            contributed = contribution_counts.get(idx, 0)
            pct_output = 100 * contributed / output_len if output_len else 0
            pct_resolved = 100 * contributed / total_filled if total_filled else 0
            data = data_counts.get(idx, 0)
            missing = missing_counts.get(idx, 0)
            ambig = ambiguity_counts.get(idx, 0)

            f.write(
                f"{idx}\t{name}\t{contributed}\t{pct_output:.2f}\t{pct_resolved:.2f}\t"
                f"{data}\t{missing}\t{ambig}\n"
            )

        f.write("\nInterpretation notes\n")
        f.write("--------------------\n")
        f.write("- The first sequence is treated as the coordinate anchor.\n")
        f.write("- Columns where sequence 1 has '-' are skipped entirely.\n")
        f.write("- This prevents lower-priority sequences from adding inserted bases caused by relaxed / incorrect mapping.\n")
        f.write("- Lower-priority sequences are only allowed to fill positions that exist in the anchor sequence.\n")
        f.write("- 'contributed_positions' = positions where this sequence was the first non-missing sequence among retained columns.\n")
        f.write("- '%_of_output' = share of the final consensus length contributed by this sequence.\n")
        f.write("- '%_of_resolved_consensus' = share among positions resolved by at least one sequence.\n")
        if treat_iupac_as_missing:
            f.write("- IUPAC ambiguity codes were treated as missing.\n")
        else:
            f.write("- IUPAC ambiguity codes were treated as data.\n")


def write_contributor_bedgraph(path, contributors):
    with open(path, "w") as f:
        f.write("chrom\tstart\tend\tcontributor_rank\n")
        if not contributors:
            return

        start = 0
        current = contributors[0]

        for i in range(1, len(contributors)):
            if contributors[i] != current:
                f.write(f"consensus\t{start}\t{i}\t{current}\n")
                start = i
                current = contributors[i]

        f.write(f"consensus\t{start}\t{len(contributors)}\t{current}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Build a priority consensus from an aligned FASTA using sequence order as priority."
    )
    parser.add_argument("input_fasta", help="Aligned input FASTA")
    parser.add_argument("output_prefix", help="Output prefix")
    parser.add_argument(
        "--treat-iupac-as-missing",
        action="store_true",
        help="Treat ambiguity codes (R,Y,S,W,K,M,B,D,H,V) as missing"
    )
    parser.add_argument(
        "--keep-insertions-relative-to-first",
        action="store_true",
        help="Do NOT ignore columns where sequence 1 has '-' (default is to ignore them)"
    )

    args = parser.parse_args()

    sequences = read_fasta(args.input_fasta)

    ignore_insertions_relative_to_first = not args.keep_insertions_relative_to_first

    results = build_priority_consensus(
        sequences,
        treat_iupac_as_missing=args.treat_iupac_as_missing,
        ignore_insertions_relative_to_first=ignore_insertions_relative_to_first
    )

    consensus_fasta = args.output_prefix + ".consensus.fasta"
    position_report = args.output_prefix + ".position_report.tsv"
    summary_report = args.output_prefix + ".summary.txt"
    contributor_track = args.output_prefix + ".contributors.bedgraph"

    write_fasta(consensus_fasta, "priority_consensus", results["consensus"])
    write_position_report(position_report, results["position_rows"], results["n_sequences"])
    write_summary_report(
        summary_report,
        sequences,
        results,
        treat_iupac_as_missing=args.treat_iupac_as_missing,
        ignore_insertions_relative_to_first=ignore_insertions_relative_to_first
    )
    write_contributor_bedgraph(contributor_track, results["contributors"])

    print(f"Read {len(sequences)} aligned sequences")
    print(f"Alignment length: {results['alignment_length']}")
    print(f"Skipped insertion columns relative to seq1: {results['skipped_insertion_columns']}")
    print(f"Output consensus length: {results['output_length']}")
    print(f"Wrote consensus FASTA: {consensus_fasta}")
    print(f"Wrote per-position report: {position_report}")
    print(f"Wrote summary report: {summary_report}")
    print(f"Wrote contributor track: {contributor_track}")


if __name__ == "__main__":
    main()
