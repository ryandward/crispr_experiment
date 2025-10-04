import argparse
from datetime import datetime
from multiprocessing import cpu_count

import gzip
import os
import json
import pandas as pd
import platform
import pysam
import re
import rich
import rich.table
import subprocess
import sys
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation
from Bio.SeqRecord import SeqRecord
from rich.console import Console
from rich.table import Table


def open_file(file, mode):
    return gzip.open(file, mode) if file.endswith(".gz") else open(file, mode)


def create_working_directory(dir_name="working_directory"):
    os.makedirs(dir_name, exist_ok=True)
    return dir_name


def create_topological_fasta(
    genbank_file_name, topological_fasta_file_name, overhang_length=0
):
    topological_records = []
    open_func = gzip.open if genbank_file_name.endswith(".gz") else open

    with open_func(genbank_file_name, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            if record.annotations.get("topology", None) == "circular":
                overhang_length = 100_000

            new_seq = record.seq + record.seq[:overhang_length]
            topological_record = SeqRecord(
                Seq(str(new_seq)),
                id=record.id,
                description=record.description,
                name=genbank_file_name,
            )
            topological_records.append(topological_record)

    with open(topological_fasta_file_name, "w") as output_handle:
        SeqIO.write(topological_records, output_handle, "fasta")


def create_fake_topological_fastq(topological_fasta_file_name, fastq_file):
    """Create a FASTQ file with phred scores from original design."""
    quality_string = (
        "I4!=======44444++++++++"  # ASCII quality scores matching original phred scores
    )

    with open(topological_fasta_file_name, "rt") as input_handle, open(
        fastq_file, "w"
    ) as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            if len(record.seq) == len(quality_string):
                record.letter_annotations["phred_quality"] = [
                    ord(c) - 33 for c in quality_string
                ]
                SeqIO.write(record, output_handle, "fastq")


def create_locus_map(genbank_file_name):
    primary_map, overhang_continue, organisms, seq_lens, topologies, all_genes = (
        {},
        {},
        {},
        {},
        {},
        {},
    )
    feature_types = {}  # Add dictionary to store feature types

    with open(genbank_file_name, "rt") as input_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            overhang_continue[record.id] = 0
            organisms[record.id] = record.annotations.get("organism", None)
            seq_lens[record.id] = len(record.seq)
            topologies[record.id] = record.annotations.get("topology", None)

            gene_count = 0
            overhang_length = 100_000 if topologies[record.id] == "circular" else 0

            for feature in record.features:
                # Get the first available qualifier as primary_id
                primary_id = None
                for key in feature.qualifiers:
                    if (
                        isinstance(feature.qualifiers[key], list)
                        and feature.qualifiers[key]
                    ):
                        primary_id = feature.qualifiers[key][0]
                        break

                if primary_id and feature.type != "gene":  # Skip 'gene' type
                    if primary_id not in feature_types:
                        feature_types[primary_id] = set()
                    feature_types[primary_id].add(feature.type)

                if feature.type == "gene":
                    gene_count += 1
                    gene_name = feature.qualifiers.get("gene", [None])[0]

                    if isinstance(feature.location, CompoundLocation):
                        # Handle wrapping genes by sorting parts by position
                        parts = sorted(feature.location.parts, key=lambda x: x.start)

                        if parts[0].start == 0 or parts[-1].end == len(record.seq):
                            # For wrapped genes, create a continuous coordinate space
                            end_segment = parts[-1]
                            start_segment = parts[0]

                            # Calculate total gene length
                            gene_length = (end_segment.end - end_segment.start) + (
                                start_segment.end - start_segment.start
                            )

                            # Map positions in first segment (after origin)
                            for position in range(
                                int(start_segment.start), int(start_segment.end)
                            ):
                                key = (record.id, position)
                                # Store the continuous coordinates
                                primary_map.setdefault(key, []).append(
                                    (
                                        primary_id,
                                        gene_name,
                                        int(
                                            end_segment.start
                                        ),  # Real start is at the end segment
                                        int(end_segment.start)
                                        + gene_length,  # Total length from real start
                                        feature.location.strand,
                                    )
                                )

                            # Map positions in second segment (before origin)
                            for position in range(
                                int(end_segment.start), int(end_segment.end)
                            ):
                                key = (record.id, position)
                                primary_map.setdefault(key, []).append(
                                    (
                                        primary_id,
                                        gene_name,
                                        int(end_segment.start),  # Same start point
                                        int(end_segment.start)
                                        + gene_length,  # Same end point
                                        feature.location.strand,
                                    )
                                )

                        else:
                            # Handle normal compound locations (not wrapping around origin)
                            for part in parts:
                                for position in range(int(part.start), int(part.end)):
                                    key = (record.id, position)
                                    primary_map.setdefault(key, []).append(
                                        (
                                            primary_id,
                                            gene_name,
                                            int(part.start),
                                            int(part.end),
                                            feature.location.strand,
                                        )
                                    )
                    else:
                        # Handle normal locations
                        for position in range(
                            int(feature.location.start), int(feature.location.end)
                        ):
                            key = (record.id, position)
                            primary_map.setdefault(key, []).append(
                                (
                                    primary_id,
                                    gene_name,
                                    int(feature.location.start),
                                    int(feature.location.end),
                                    feature.location.strand,
                                )
                            )

                            # Add overhang positions for circular genomes
                            if (
                                topologies[record.id] == "circular"
                                and position < overhang_length
                            ):
                                key = (record.id, position + len(record.seq))
                                primary_map.setdefault(key, []).append(
                                    (
                                        primary_id,
                                        gene_name,
                                        int(feature.location.start) + len(record.seq),
                                        int(feature.location.end) + len(record.seq),
                                        feature.location.strand,
                                    )
                                )

            all_genes[record.id] = gene_count

    return primary_map, organisms, seq_lens, topologies, all_genes, feature_types


def get_true_chrom_lengths(gb_file):
    chrom_lengths = {}
    with open(gb_file, "r") as f:
        for rec in SeqIO.parse(f, "genbank"):
            chrom_lengths[rec.id] = len(rec.seq)
    return chrom_lengths


def get_topological_chrom_lengths(fasta_file):
    chrom_lengths = {}
    with open(fasta_file, "r") as f:
        for rec in SeqIO.parse(f, "fasta"):
            chrom_lengths[rec.id] = len(rec.seq)
    return chrom_lengths


def get_diff(spacer, target):
    differences = [
        f"{target_nt}{i + 1}{spacer_nt}"
        for i, (target_nt, spacer_nt) in enumerate(zip(target, spacer))
        if target_nt != spacer_nt
    ]
    return ",".join(differences) if differences else None


def get_coords(tar_start, tar_end, chrom_length):
    start_circular = tar_start % chrom_length
    end_circular = (
        tar_end % chrom_length if tar_end % chrom_length != 0 else chrom_length
    )
    return (
        f"({start_circular}..{chrom_length}, 0..{end_circular})"
        if start_circular > end_circular
        else f"{start_circular}..{end_circular}"
    )


def get_offset(
    target_dir, tar_start, tar_end, feature_start, feature_end, chrom_length
):
    # First normalize all coordinates to be within genome length
    tar_start = tar_start % chrom_length
    tar_end = tar_end % chrom_length
    feature_start = feature_start % chrom_length
    feature_end = feature_end % chrom_length

    if target_dir == "F":
        if feature_start > feature_end:  # Gene wraps around origin
            if tar_start < feature_end:  # Target is in early part
                return tar_start + (chrom_length - feature_start)
            else:  # Target is in later part
                return tar_start - feature_start
        else:  # Normal gene
            return tar_start - feature_start
    if target_dir == "R":
        if feature_start > feature_end:  # Gene wraps around origin
            if tar_end < feature_end:  # Target is in early part
                return (chrom_length - tar_end) + (feature_end - feature_start)
            else:  # Target is in later part
                return feature_end - tar_end
        else:  # Normal gene
            return feature_end - tar_end
    return None


def get_overlap(tar_start, tar_end, feature_start, feature_end, chrom_length):
    # First normalize all coordinates to be within genome length
    tar_start = tar_start % chrom_length
    tar_end = tar_end % chrom_length
    feature_start = feature_start % chrom_length
    feature_end = feature_end % chrom_length

    # Handle cases where coordinates cross the origin
    if feature_start > feature_end:  # Feature wraps around origin
        if tar_start < feature_end:  # Target is in the start region
            overlap_start = tar_start
            overlap_end = min(tar_end, feature_end)
        elif tar_start >= feature_start:  # Target is in the end region
            overlap_start = tar_start
            overlap_end = min(tar_end, chrom_length)
            if tar_end < feature_end:  # Add additional overlap in start region
                overlap_end += min(tar_end, feature_end)
        else:
            return 0
    else:  # Normal feature
        overlap_start = max(tar_start, feature_start)
        overlap_end = min(tar_end, feature_end)

    overlap = overlap_end - overlap_start
    return overlap if overlap >= 0 else 0


def convert_degenerate_to_regex(pam_pattern):
    """
    Convert IUPAC degenerate nucleotide codes to regex character classes.
    """
    degenerate_map = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "R": "[AG]",  # puRines
        "Y": "[CT]",  # pYrimidines
        "S": "[CG]",  # Strong (3 H-bonds)
        "W": "[AT]",  # Weak (2 H-bonds)
        "K": "[GT]",  # Keto
        "M": "[AC]",  # aMino
        "B": "[CGT]",  # not A
        "D": "[AGT]",  # not C
        "H": "[ACT]",  # not G
        "V": "[ACG]",  # not T
        "N": ".",  # any nucleotide
    }

    regex_pattern = ""
    for char in pam_pattern.upper():
        regex_pattern += degenerate_map.get(char, char)

    return regex_pattern


def pam_matches(pam_pattern, extracted_pam):
    if not extracted_pam:
        return False
    if pam_pattern == "N" * len(pam_pattern) or not pam_pattern:
        return True

    # Only check the first len(pam_pattern) characters of extracted_pam
    if len(extracted_pam) < len(pam_pattern):
        return False

    test_seq = extracted_pam[: len(pam_pattern)]
    regex_pattern = convert_degenerate_to_regex(pam_pattern)
    return bool(re.match(regex_pattern + "$", test_seq))


def process_single_alignment(
    read,
    primary_map,
    fasta,
    true_chrom_lengths,
    pam,
    pam_direction,
    regulatory,
):
    """Process a single alignment and return a dict of results"""
    if read.query_sequence is None:
        return None

    # For downstream PAM (e.g., NGG):
    # Forward reads already have CCG at start
    # Reverse reads need to be rev-comped to get CCG at start
    # if pam_direction == "downstream":
    if read.is_reverse:
        full_sequence = str(Seq(read.query_sequence))
    else:
        full_sequence = str(Seq(read.query_sequence).reverse_complement())
    # else:
    #     full_sequence = read.query_sequence

    # Extract PAM and spacer
    extracted_pam = (
        full_sequence[-len(pam) :] if len(full_sequence) >= len(pam) else None
    )
    spacer = (
        full_sequence[: -len(pam)] if len(full_sequence) > len(pam) else full_sequence
    )

    rows = {
        "name": read.query_name,
        "spacer": spacer,
        "len": len(spacer),
    }

    if read.is_unmapped:
        return rows

    if read.is_mapped:
        try:
            if read.is_reverse:
                target = str(Seq(read.get_reference_sequence()))
            else:
                target = str(Seq(read.get_reference_sequence()).reverse_complement())

            target = target[: -len(pam)]  # Trim PAM

            chr = read.reference_name
            chrom_length = true_chrom_lengths.get(chr, None)

            # Coordinate handling stays the same
            if pam_direction == "downstream":
                if not read.is_reverse:
                    tar_start = (read.reference_start + len(pam)) % chrom_length
                    tar_end = read.reference_end % chrom_length
                else:
                    tar_start = read.reference_start % chrom_length
                    tar_end = (read.reference_end - len(pam)) % chrom_length
            else:
                tar_start = read.reference_start % chrom_length
                tar_end = (read.reference_end - len(pam)) % chrom_length

            sp_dir = "F" if not read.is_reverse else "R"
            coords = get_coords(tar_start, tar_end, chrom_length)

            type = "mismatch" if target != rows["spacer"] else "perfect"
            diff = get_diff(rows["spacer"], target)

        except ValueError as e:
            console.log(f"Error processing read {read.query_name}: {e}")
            sys.exit(1)

        rows.update(
            {
                "target": target,
                "chr": chr,
                "tar_start": tar_start,
                "tar_end": tar_end,
                "sp_dir": sp_dir,
                "pam": extracted_pam,
                "coords": coords,
                "type": type,
                "diff": diff,
            }
        )

        if sp_dir == "F":
            aligned_genes = {
                gene_info
                for pos in range(tar_start + regulatory, tar_end + regulatory)
                for gene_info in primary_map.get((chr, pos), [])
            }
        else:
            aligned_genes = {
                gene_info
                for pos in range(tar_start - regulatory, tar_end - regulatory)
                for gene_info in primary_map.get((chr, pos), [])
            }

        if not aligned_genes:
            rows.update(
                {
                    "primary_id": None,
                    "gene": None,
                    "offset": None,
                    "overlap": None,
                    "tar_dir": None,
                }
            )
            return rows
        else:
            for (
                primary_id,
                gene_name,
                feature_start,
                feature_end,
                feature_strand,
            ) in aligned_genes:
                target_orientation = (
                    "F"
                    if feature_strand == 1
                    else "R" if feature_strand == -1 else None
                )
                offset = get_offset(
                    target_orientation,
                    tar_start,
                    tar_end,
                    feature_start,
                    feature_end,
                    chrom_length,
                )
                overlap = get_overlap(
                    tar_start,
                    tar_end,
                    feature_start,
                    feature_end,
                    chrom_length,
                )

                rows_copy = rows.copy()
                rows_copy.update(
                    {
                        "primary_id": primary_id,
                        "gene": gene_name if gene_name else primary_id,
                        "offset": offset,
                        "overlap": overlap,
                        "tar_dir": target_orientation,
                    }
                )
                return rows_copy


def parse_sam_output_streaming(
    samfile,
    primary_map,
    topological_fasta_file_name,
    gb_file_name,
    pam,
    pam_direction,
    regulatory,
):
    """Stream process alignments one at a time"""
    true_chrom_lengths = get_true_chrom_lengths(gb_file_name)
    console = Console(file=sys.stderr)
    processed = 0

    # Write header to stdout, making sure to flush
    header = [
        "spacer",
        "primary_id",
        "gene",
        "feature_types",
        "chr",
        "pam",
        "target",
        "tar_start",
        "tar_end",
        "offset",
        "overlap",
        "sp_dir",
        "tar_dir",
    ]
    sys.stdout.write("\t".join(header) + "\n")
    sys.stdout.flush()

    with pysam.FastaFile(topological_fasta_file_name) as fasta:
        for read in samfile.fetch(until_eof=True):
            processed += 1
            if processed % 10000 == 0:
                console.log(f"Processed {processed:,} alignments")
                sys.stderr.flush()

            result = process_single_alignment(
                read,
                primary_map,
                fasta,
                true_chrom_lengths,
                pam,
                pam_direction,
                regulatory,
            )

            if result:
                # Write directly to stdout and flush
                line = "\t".join(str(result.get(col, "None")) for col in header)
                sys.stdout.write(line + "\n")
                sys.stdout.flush()


def filter_offtargets_by_pam(df):
    targeting_spacers = df[df["target"].notna()]["spacer"].unique()
    return df[~((df["target"].isna()) & (df["spacer"].isin(targeting_spacers)))]


def main(args):
    console = Console(file=sys.stderr)
    console.log("[bold red]Initializing barcode target seeker[/bold red]")

    # Use a temp directory with more space than /tmp
    temp_base = (
        os.path.expanduser("~/tmp")
        if os.path.exists(os.path.expanduser("~/tmp"))
        else None
    )
    with tempfile.TemporaryDirectory(dir=temp_base) as working_dir:
        # First, handle BAM conversion and sorting
        console.log("Converting SAM to sorted BAM...")
        bam_file_path = os.path.join(working_dir, "alignment.bam")
        sorted_bam_file_path = os.path.join(working_dir, "alignment.sorted.bam")

        subprocess.run(
            ["samtools", "view", "-bS", args.sam_file, "-o", bam_file_path], check=True
        )
        subprocess.run(
            ["samtools", "sort", bam_file_path, "-o", sorted_bam_file_path], check=True
        )
        subprocess.run(["samtools", "index", sorted_bam_file_path], check=True)

        # Now create topological FASTA
        console.log("Creating topological FASTA reference...")
        topological_fasta_file_name = os.path.join(
            working_dir,
            os.path.splitext(os.path.basename(args.genome_file))[0] + ".fasta",
        )
        create_topological_fasta(args.genome_file, topological_fasta_file_name)

        # Generate locus map
        console.log("Generating gene feature maps...")
        primary_map, organisms, seq_lens, topologies, all_genes, feature_types = (
            create_locus_map(args.genome_file)
        )

        # Stream process alignments
        console.log("Processing alignments...")
        with pysam.AlignmentFile(sorted_bam_file_path, "rb") as samfile:
            parse_sam_output_streaming(
                samfile,
                primary_map,
                topological_fasta_file_name,
                args.genome_file,
                args.pam,
                args.pam_direction,
                args.regulatory,
            )

        console.log("[bold green]Processing complete[/bold green]", file=sys.stderr)
        sys.stderr.flush()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Map barcodes to a circular genome",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("sam_file", help="Path to an existing SAM file", type=str)
    parser.add_argument("genome_file", help="Path to genome_gb_file", type=str)
    parser.add_argument("pam", help="PAM sequence", type=str)
    parser.add_argument(
        "--pam-direction",
        choices=["upstream", "downstream"],
        default="downstream",
        help="Direction of the PAM sequence",
    )
    parser.add_argument(
        "--regulatory",
        type=int,
        default=0,
        help="Regulatory region size around target position.",
    )
    args = parser.parse_args()
    main(args)
