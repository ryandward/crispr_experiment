import argparse
from datetime import datetime
from multiprocessing import cpu_count
import gc  # Add garbage collection

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
    if topological_fasta_file_name.endswith(".gz"):
        with gzip.open(topological_fasta_file_name, "rt") as input_handle, open(
            fastq_file, "w"
        ) as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")
    else:
        with open(topological_fasta_file_name, "rt") as input_handle, open(
            fastq_file, "w"
        ) as output_handle:
            for record in SeqIO.parse(input_handle, "fasta"):
                record.letter_annotations["phred_quality"] = [40] * len(record)
                SeqIO.write(record, output_handle, "fastq")


def create_locus_map(genbank_file_name):
    locus_map, overhang_continue, organisms, seq_lens, topologies, all_genes = (
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
                locus_tag = feature.qualifiers.get("locus_tag", [None])[0]
                if locus_tag and feature.type != "gene":  # Skip 'gene' type
                    if locus_tag not in feature_types:
                        feature_types[locus_tag] = set()
                    feature_types[locus_tag].add(feature.type)

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
                                locus_map.setdefault(key, []).append(
                                    (
                                        locus_tag,
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
                                locus_map.setdefault(key, []).append(
                                    (
                                        locus_tag,
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
                                    locus_map.setdefault(key, []).append(
                                        (
                                            locus_tag,
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
                            locus_map.setdefault(key, []).append(
                                (
                                    locus_tag,
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
                                locus_map.setdefault(key, []).append(
                                    (
                                        locus_tag,
                                        gene_name,
                                        int(feature.location.start) + len(record.seq),
                                        int(feature.location.end) + len(record.seq),
                                        feature.location.strand,
                                    )
                                )

            all_genes[record.id] = gene_count

    return locus_map, organisms, seq_lens, topologies, all_genes, feature_types


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

    IUPAC codes:
    A = Adenine
    C = Cytosine
    G = Guanine
    T = Thymine
    R = A or G (puRines)
    Y = C or T (pYrimidines)
    S = C or G (Strong)
    W = A or T (Weak)
    K = G or T (Keto)
    M = A or C (aMino)
    B = C, G, or T (not A)
    D = A, G, or T (not C)
    H = A, C, or T (not G)
    V = A, C, or G (not T)
    N = any nucleotide
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


def extract_downstream_pam(
    pam,
    tar_start,
    tar_end,
    chrom,
    fasta,
    dir,
    true_chrom_lengths,
    topological_chrom_lengths,
):
    true_chrom_length = true_chrom_lengths.get(chrom, None)
    topological_chrom_length = topological_chrom_lengths.get(chrom, None)

    if not pam or None in (
        tar_start,
        tar_end,
        chrom,
        fasta,
        dir,
        true_chrom_length,
        topological_chrom_length,
    ):
        return None

    if dir == "F":
        if tar_end + len(pam) > topological_chrom_length:
            # Handle wrapping PAM by concatenating end and start of sequence
            end_part = fasta.fetch(reference=chrom, start=tar_end).upper()
            start_part = fasta.fetch(
                reference=chrom, start=0, end=(tar_end + len(pam)) % true_chrom_length
            ).upper()
            extracted_pam = end_part + start_part
            return extracted_pam[: len(pam)]
        extracted_pam = fasta.fetch(
            reference=chrom, start=tar_end, end=tar_end + len(pam)
        ).upper()
    elif dir == "R":
        if tar_start - len(pam) < 0:
            return None
        extracted_pam = fasta.fetch(
            reference=chrom, start=tar_start - len(pam), end=tar_start
        ).upper()
        extracted_pam = str(Seq(extracted_pam).reverse_complement())
    else:
        return None
    return extracted_pam


def extract_upstream_pam(
    pam,
    tar_start,
    tar_end,
    chrom,
    fasta,
    dir,
    true_chrom_lengths,
    topological_chrom_lengths,
):
    true_chrom_length = true_chrom_lengths.get(chrom, None)
    topological_chrom_length = topological_chrom_lengths.get(chrom, None)
    if not pam or None in (
        tar_start,
        tar_end,
        chrom,
        fasta,
        dir,
        true_chrom_length,
        topological_chrom_length,
    ):
        return None

    if dir == "F":
        if tar_start - len(pam) < 0:
            # Handle wrapping PAM by concatenating end and start of sequence
            end_part = fasta.fetch(
                reference=chrom, start=true_chrom_length + (tar_start - len(pam))
            ).upper()
            start_part = fasta.fetch(reference=chrom, start=0, end=tar_start).upper()
            extracted_pam = end_part + start_part
            return extracted_pam[-len(pam) :]
        extracted_pam = fasta.fetch(
            reference=chrom, start=tar_start - len(pam), end=tar_start
        ).upper()
    elif dir == "R":
        if tar_end + len(pam) > topological_chrom_length:
            return None
        extracted_pam = fasta.fetch(
            reference=chrom, start=tar_end, end=tar_end + len(pam)
        ).upper()
        extracted_pam = str(Seq(extracted_pam).reverse_complement())
    else:
        return None
    return extracted_pam


def parse_sam_output(
    samfile,
    locus_map,
    topological_fasta_file_name,
    gb_file_name,
    pam,
    pam_direction,
    regulatory,
):
    rows_list = []
    true_chrom_lengths = get_true_chrom_lengths(gb_file_name)
    topological_chrom_lengths = get_topological_chrom_lengths(
        topological_fasta_file_name
    )

    with pysam.FastaFile(topological_fasta_file_name) as fasta:
        for read in samfile.fetch():
            if read.query_sequence is None:
                continue

            extracted_pam = None

            if pam:
                if pam_direction == "downstream":
                    extracted_pam = extract_downstream_pam(
                        pam,
                        read.reference_start,
                        read.reference_end,
                        read.reference_name,
                        fasta,
                        "F" if not read.is_reverse else "R",
                        true_chrom_lengths,
                        topological_chrom_lengths,
                    )
                elif pam_direction == "upstream":
                    extracted_pam = extract_upstream_pam(
                        pam,
                        read.reference_start,
                        read.reference_end,
                        read.reference_name,
                        fasta,
                        "F" if not read.is_reverse else "R",
                        true_chrom_lengths,
                        topological_chrom_lengths,
                    )

                if read.is_mapped and not pam_matches(pam, extracted_pam):
                    read.is_unmapped = True
                    extracted_pam = None

            rows = {
                "name": read.query_name,
                "spacer": (
                    read.query_sequence
                    if not read.is_reverse
                    else str(Seq(read.query_sequence).reverse_complement())
                ),
                "len": len(read.query_sequence),
            }

            if read.is_unmapped:
                rows_list.append(rows)
                continue

            if read.is_mapped:
                try:
                    target = (
                        read.get_reference_sequence()
                        if not read.is_reverse
                        else str(
                            Seq(read.get_reference_sequence()).reverse_complement()
                        )
                    )
                    mismatches = int(read.get_tag("NM"))
                    chr = read.reference_name
                    chrom_length = true_chrom_lengths.get(chr, None)

                    # Always use modulo for consistent positive coordinates
                    tar_start = read.reference_start % chrom_length
                    tar_end = read.reference_end % chrom_length

                    sp_dir = "F" if not read.is_reverse else "R"
                    coords = get_coords(tar_start, tar_end, chrom_length)
                    type = "mismatch" if mismatches > 0 else "perfect"
                    diff = get_diff(rows["spacer"], target)

                except ValueError as e:
                    print(e, file=sys.stderr)
                    sys.exit(1)

                rows.update(
                    {
                        "target": target,
                        "mismatches": mismatches,
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
                        for gene_info in locus_map.get((chr, pos), [])
                    }
                else:
                    aligned_genes = {
                        gene_info
                        for pos in range(tar_start - regulatory, tar_end - regulatory)
                        for gene_info in locus_map.get((chr, pos), [])
                    }

                if not aligned_genes:
                    rows.update(
                        {
                            "locus_tag": None,
                            "gene": None,
                            "offset": None,
                            "overlap": None,
                            "tar_dir": None,
                        }
                    )
                    rows_list.append(rows)
                else:
                    for (
                        locus_tag,
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
                            chrom_length,  # Pass chrom_length to get_offset
                        )
                        overlap = get_overlap(
                            tar_start, tar_end, feature_start, feature_end, chrom_length
                        )

                        rows_copy = rows.copy()
                        rows_copy.update(
                            {
                                "locus_tag": locus_tag,
                                "gene": gene_name if gene_name else locus_tag,
                                "offset": offset,
                                "overlap": overlap,
                                "tar_dir": target_orientation,
                            }
                        )
                        rows_list.append(rows_copy)

    return rows_list


def run_bowtie_and_parse(
    sgrna_fastq_file_name,
    topological_fasta_file_name,
    locus_map,
    num_mismatches,
    num_threads,
    pam,
    pam_direction,
    regulatory,
):
    results = []
    # Use a temp directory with more space than /tmp
    temp_base = (
        os.path.expanduser("~/tmp")
        if os.path.exists(os.path.expanduser("~/tmp"))
        else None
    )
    with tempfile.TemporaryDirectory(dir=temp_base) as temp_dir:
        # Create filenames with shorter paths to avoid C++ string length issues
        index_prefix = os.path.join(temp_dir, "genome")

        # Fix the typo in command name
        print("Building bowtie index...", file=sys.stderr)
        bowtie_build_command = [
            "bowtie-build",  # This was incorrectly "bowtie-buAre ild"
            "--quiet",
            topological_fasta_file_name,
            index_prefix,
        ]

        # Use subprocess.run for better error handling
        try:
            subprocess.run(bowtie_build_command, check=True, stderr=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            print(f"Bowtie-build failed: {e.stderr}", file=sys.stderr)
            raise RuntimeError("Failed to build bowtie index")

        # Write SAM directly to file to avoid C++ stream issues
        sam_file = os.path.join(temp_dir, "alignment.sam")

        # Run bowtie with explicit output file
        bowtie_command = [
            "bowtie",
            "-S",
            "-k",
            "1",  # Only unique mappings - no off-targets
            "--best",
            "-v",
            str(num_mismatches),
            "-p",
            str(num_threads),
            index_prefix,
            sgrna_fastq_file_name,
            sam_file,  # Write to file instead of pipe
        ]

        print(f"Running: {' '.join(bowtie_command)}", file=sys.stderr)

        try:
            process = subprocess.run(
                bowtie_command, check=False, stderr=subprocess.PIPE, text=True
            )
            if process.returncode != 0:
                print(f"Bowtie stderr: {process.stderr}", file=sys.stderr)
                raise RuntimeError("Bowtie alignment failed")
        except Exception as e:
            print(f"Bowtie execution error: {str(e)}", file=sys.stderr)
            raise

        # Now read from the SAM file
        try:
            with pysam.AlignmentFile(sam_file, "r") as samfile:
                results = parse_sam_output(
                    samfile,
                    locus_map,
                    topological_fasta_file_name,
                    args.genome_file,
                    pam,
                    pam_direction,
                    regulatory,
                )
        except Exception as e:
            print(f"Error processing SAM file: {str(e)}", file=sys.stderr)
            raise

    if not results:
        raise RuntimeError(
            "No results were returned from the Bowtie process. Check your input files and parameters."
        )
    return results


def filter_offtargets_by_pam(df):
    try:
        import polars as pl

        # Get spacers that have valid targets (Polars syntax)
        targeting_spacers = (
            df.filter(pl.col("target").is_not_null())["spacer"].unique().to_list()
        )

        # Filter out rows where target is null AND spacer is in targeting_spacers list
        return df.filter(
            ~(pl.col("target").is_null() & pl.col("spacer").is_in(targeting_spacers))
        )
    except ImportError:
        # Pandas fallback
        targeting_spacers = df[df["target"].notna()]["spacer"].unique().tolist()
        return df[~(df["target"].isna() & df["spacer"].isin(targeting_spacers))]


def main(args):
    # Echo all received arguments to stderr for debugging
    sys.stderr.write(f"Received arguments: {args}\n")
    sys.stderr.flush()

    console = Console(file=sys.stderr)
    console.log("[bold red]Initializing barcode target seeker[/bold red]")
    num_threads = cpu_count() // 2
    # Use a temp directory with more space than /tmp
    temp_base = (
        os.path.expanduser("~/tmp")
        if os.path.exists(os.path.expanduser("~/tmp"))
        else None
    )
    with tempfile.TemporaryDirectory(dir=temp_base) as working_dir:
        topological_fasta_file_name = os.path.join(
            working_dir,
            os.path.splitext(os.path.basename(args.genome_file))[0] + ".fasta",
        )

        console.log("Annotating regions to identify...")

        base_name = os.path.basename(args.sgrna_file)
        if ".fastq" in base_name:
            sgrna_fastq_file_name = args.sgrna_file
        elif base_name.endswith(".fasta") or base_name.endswith(".fasta.gz"):
            ext = ".fasta.gz" if base_name.endswith(".fasta.gz") else ".fasta"
            base_name = base_name[: -len(ext)]
            sgrna_fastq_file_name = os.path.join(working_dir, base_name + ".fastq")
            create_fake_topological_fastq(args.sgrna_file, sgrna_fastq_file_name)
        else:
            console.log(f"[bold red]File extension not recognized. [/bold red]")
            sys.exit(1)

        console.log("Generating topological coordinate maps...")
        locus_map, organisms, seq_lens, topologies, all_genes, feature_types = (
            create_locus_map(args.genome_file)
        )

        create_topological_fasta(args.genome_file, topological_fasta_file_name)

        console.log("Aligning annotations to genome...")
        results = run_bowtie_and_parse(
            sgrna_fastq_file_name,
            topological_fasta_file_name,
            locus_map,
            args.mismatches,
            num_threads,
            args.pam,
            args.pam_direction,
            args.regulatory,
        )

        with pysam.FastaFile(topological_fasta_file_name) as _:
            pass

    try:
        console.log("Finding matches...")
        console.log(f"Processing {len(results):,} raw alignment results...")

        # Stream processing without loading all results into memory
        console.log("Starting streaming deduplication...")

        try:
            import polars as pl

            # Stream process in very small batches to avoid memory issues
            batch_size = 25_000  # Reduced batch size for better memory management
            processed_count = 0
            output_file = tempfile.NamedTemporaryFile(
                mode="w", delete=False, suffix=".tsv"
            )
            header_written = False

            # Use a much smaller memory-efficient deduplication approach
            # Only keep track of recent duplicates, not all seen keys
            recent_keys = {}
            max_recent_keys = 100_000  # Limit memory usage

            try:
                for i in range(0, len(results), batch_size):
                    batch = results[i : i + batch_size]
                    processed_count += len(batch)

                    if processed_count % 1_000_000 == 0:
                        console.log(
                            f"Streaming processed {processed_count:,} results..."
                        )

                    # Quick deduplication within batch only
                    batch_keys_seen = set()
                    batch_dedupe = []
                    for result in batch:
                        key = (
                            result.get("chr"),
                            result.get("tar_start"),
                            result.get("tar_end"),
                            result.get("spacer"),
                            result.get("target"),
                        )
                        # Only deduplicate within this batch and recent batches
                        if key not in batch_keys_seen and key not in recent_keys:
                            batch_keys_seen.add(key)
                            batch_dedupe.append(result)

                    # Update recent keys with limited memory
                    for key in batch_keys_seen:
                        recent_keys[key] = processed_count

                    # Periodically clean old keys to prevent memory growth
                    if len(recent_keys) > max_recent_keys:
                        # Remove oldest 20% of keys
                        sorted_keys = sorted(recent_keys.items(), key=lambda x: x[1])
                        keys_to_remove = sorted_keys[: len(sorted_keys) // 5]
                        for key, _ in keys_to_remove:
                            del recent_keys[key]

                    if not batch_dedupe:
                        continue

                    # Process this small batch with polars
                    # Define explicit schema to avoid type inference issues
                    schema = {
                        "name": pl.Utf8,
                        "spacer": pl.Utf8,
                        "len": pl.Int64,
                        "target": pl.Utf8,
                        "mismatches": pl.Int64,
                        "chr": pl.Utf8,
                        "tar_start": pl.Int64,
                        "tar_end": pl.Int64,
                        "sp_dir": pl.Utf8,
                        "pam": pl.Utf8,
                        "coords": pl.Utf8,
                        "type": pl.Utf8,
                        "diff": pl.Utf8,
                        "locus_tag": pl.Utf8,
                        "gene": pl.Utf8,
                        "offset": pl.Int64,
                        "overlap": pl.Int64,
                        "tar_dir": pl.Utf8,
                    }
                    batch_df = pl.DataFrame(batch_dedupe, schema=schema)

                    # Add feature types
                    batch_df = batch_df.with_columns(
                        [
                            pl.col("locus_tag")
                            .map_elements(
                                lambda x: (
                                    ",".join(sorted(feature_types.get(x, [])))
                                    if x
                                    else None
                                ),
                                return_dtype=pl.Utf8,
                            )
                            .alias("feature_types")
                        ]
                    )

                    # Filter offtargets by PAM
                    batch_df = filter_offtargets_by_pam(batch_df)

                    if batch_df.height == 0:
                        continue

                    # Convert to pandas for final processing
                    batch_pd = batch_df.to_pandas()

                    # Add site column
                    batch_pd.loc[batch_pd["target"].notnull(), "site"] = (
                        batch_pd["chr"].astype(str)
                        + "_"
                        + batch_pd["tar_start"].astype(str)
                        + ".."
                        + batch_pd["tar_end"].astype(str)
                    )

                    # Calculate stats for this batch only
                    site_counts = batch_pd.groupby("spacer")["site"].nunique()
                    gene_counts = batch_pd.groupby("spacer")["gene"].nunique(
                        dropna=True
                    )
                    intergenic_counts = (
                        batch_pd.loc[batch_pd["gene"].isna()]
                        .groupby("spacer")["site"]
                        .nunique()
                    )

                    batch_stats = (
                        pd.DataFrame(
                            {
                                "sites": site_counts,
                                "genes": gene_counts,
                                "intergenic": intergenic_counts,
                            }
                        )
                        .fillna(0)
                        .astype(int)
                    )

                    batch_pd = batch_pd.merge(
                        batch_stats, left_on="spacer", right_index=True, how="left"
                    )

                    # Column ordering
                    column_order = [
                        "spacer",
                        "locus_tag",
                        "gene",
                        "feature_types",
                        "chr",
                    ]
                    if not (
                        batch_pd["pam"].isnull().all() or batch_pd["pam"].nunique() == 1
                    ):
                        column_order.append("pam")
                    column_order.append("mismatches")
                    column_order.extend(
                        [
                            "target",
                            "tar_start",
                            "tar_end",
                            "offset",
                            "overlap",
                            "sp_dir",
                            "tar_dir",
                            "sites",
                            "genes",
                            "intergenic",
                        ]
                    )

                    batch_pd = batch_pd[column_order]

                    # Ensure integer columns are properly formatted (no scientific notation)
                    integer_columns = [
                        "tar_start",
                        "tar_end",
                        "offset",
                        "overlap",
                        "mismatches",
                        "sites",
                        "genes",
                        "intergenic",
                    ]
                    for col in integer_columns:
                        if col in batch_pd.columns:
                            batch_pd[col] = batch_pd[col].astype(
                                "Int64"
                            )  # Nullable integer type

                    # Write batch to temp file
                    if not header_written:
                        batch_pd.to_csv(
                            output_file,
                            sep="\t",
                            index=False,
                            na_rep="None",
                            float_format="%.0f",
                        )
                        header_written = True
                    else:
                        batch_pd.to_csv(
                            output_file,
                            sep="\t",
                            index=False,
                            header=False,
                            na_rep="None",
                            float_format="%.0f",
                        )

                    # Clear memory aggressively
                    del batch_df, batch_pd, batch_dedupe
                    # Also clear the batch data
                    batch = None

                    # Force garbage collection every few batches to prevent memory buildup
                    if processed_count % 500_000 == 0:
                        gc.collect()

                output_file.close()

                console.log(
                    f"✓ Streaming complete: {processed_count:,} results processed"
                )

                console.log("Sorting results by genomic coordinates...")

                # Sort the file by chromosome and position using system sort for efficiency
                import subprocess

                sorted_file = tempfile.NamedTemporaryFile(
                    mode="w", delete=False, suffix=".tsv"
                )
                sorted_file.close()

                try:
                    # Use system sort: -k5,5 (chromosome) -k9,9n (tar_start numeric) -k10,10n (tar_end numeric)
                    sort_cmd = [
                        "sort",
                        "-t",
                        "\t",
                        "-k5,5",  # Sort by chromosome (5th column)
                        "-k9,9n",  # Then by tar_start numerically (9th column)
                        "-k10,10n",  # Then by tar_end numerically (10th column)
                        output_file.name,
                    ]

                    with open(sorted_file.name, "w") as f:
                        subprocess.run(sort_cmd, stdout=f, check=True)

                    # Stream the sorted file to stdout
                    with open(sorted_file.name, "r") as f:
                        for line in f:
                            sys.stdout.write(line)

                    # Clean up temp files
                    os.unlink(sorted_file.name)

                except Exception as sort_error:
                    console.log(
                        f"[yellow]Warning: Could not sort results: {sort_error}[/yellow]"
                    )
                    console.log("Outputting unsorted results...")
                    # Fallback to unsorted output
                    with open(output_file.name, "r") as f:
                        for line in f:
                            sys.stdout.write(line)

                # Clean up temp file
                os.unlink(output_file.name)

                # Simple stats for summary table
                results = pd.DataFrame(
                    {
                        "spacer": ["dummy"],
                        "chr": ["dummy"],
                        "locus_tag": [None],
                        "target": [None],
                    }
                )
                off_target_spacers = 0

            except Exception as e:
                console.log(f"[bold red]Error during streaming: {e}[/bold red]")
                if (
                    "output_file" in locals()
                    and hasattr(output_file, "name")
                    and os.path.exists(output_file.name)
                ):
                    os.unlink(output_file.name)
                raise

        except ImportError:
            console.log(
                "[bold red]Polars not available, using pandas streaming[/bold red]"
            )
            # Pandas fallback with streaming
            processed_count = 0
            seen_keys = set()
            output_chunks = []

            for i in range(0, len(results), 25_000):  # Smaller batches for pandas
                batch = results[i : i + 25_000]
                processed_count += len(batch)

                if processed_count % 500_000 == 0:
                    console.log(f"Pandas streaming {processed_count:,} results...")

                # Quick dedup
                batch_dedupe = []
                for result in batch:
                    key = (
                        result.get("chr"),
                        result.get("tar_start"),
                        result.get("tar_end"),
                        result.get("spacer"),
                        result.get("target"),
                    )
                    if key not in seen_keys:
                        seen_keys.add(key)
                        batch_dedupe.append(result)

                if batch_dedupe:
                    output_chunks.extend(batch_dedupe)

                # Process accumulated chunks when they get too big
                if len(output_chunks) > 100_000:
                    chunk_df = pd.DataFrame(output_chunks)
                    # Minimal processing - just output what we have
                    column_order = [
                        "spacer",
                        "locus_tag",
                        "gene",
                        "chr",
                        "pam",
                        "mismatches",
                        "target",
                        "tar_start",
                        "tar_end",
                        "offset",
                        "overlap",
                        "sp_dir",
                        "tar_dir",
                    ]
                    available_cols = [
                        col for col in column_order if col in chunk_df.columns
                    ]

                    # Ensure integer columns are properly formatted
                    integer_columns = [
                        "tar_start",
                        "tar_end",
                        "offset",
                        "overlap",
                        "mismatches",
                    ]
                    for col in integer_columns:
                        if col in chunk_df.columns:
                            chunk_df[col] = chunk_df[col].astype("Int64")

                    chunk_df[available_cols].to_csv(
                        sys.stdout,
                        sep="\t",
                        index=False,
                        na_rep="None",
                        header=(processed_count <= 100_000),
                        float_format="%.0f",
                    )
                    output_chunks = []
                    del chunk_df

            # Output remaining chunks
            if output_chunks:
                chunk_df = pd.DataFrame(output_chunks)
                available_cols = [
                    col for col in column_order if col in chunk_df.columns
                ]

                # Ensure integer columns are properly formatted
                integer_columns = [
                    "tar_start",
                    "tar_end",
                    "offset",
                    "overlap",
                    "mismatches",
                ]
                for col in integer_columns:
                    if col in chunk_df.columns:
                        chunk_df[col] = chunk_df[col].astype("Int64")

                chunk_df[available_cols].to_csv(
                    sys.stdout,
                    sep="\t",
                    index=False,
                    na_rep="None",
                    header=(len(output_chunks) == len(results)),
                    float_format="%.0f",
                )

            console.log(
                f"✓ Pandas streaming complete: {processed_count:,} -> {len(seen_keys):,} unique results"
            )
            return  # Skip the rest since we already output data

    except FileNotFoundError:
        console.log(
            f"[bold red]Trouble with Bowtie aligner. Try using a lower number of mismatches.[/bold red]"
        )
        sys.exit(1)
    except KeyError as e:
        console.log(
            f"[bold red]All of the proposed barcodes are missing some key attributes[/bold red]: {e}"
        )
        sys.exit(1)

    combined_table = Table(
        box=rich.table.box.SIMPLE_HEAVY,
        caption=f"Analysis completed at [u]{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}[/u]",
        highlight=True,
        show_header=True,
        style="white",
    )
    combined_table.add_column("Category", justify="right", style="white", min_width=30)
    combined_table.add_column("Details", justify="left", style="white", min_width=40)

    # Input Configuration Section
    combined_table.add_section()
    combined_table.add_row(
        "[bold bright_magenta]Input Configuration[/bold bright_magenta]",
        "[italic]Analysis Parameters[/italic]",
    )
    combined_table.add_row("Barcodes File", os.path.abspath(args.sgrna_file))
    combined_table.add_row("Genome File", os.path.abspath(args.genome_file))
    combined_table.add_row("PAM Sequence", args.pam)
    combined_table.add_row(
        "PAM Direction",
        "Downstream" if args.pam_direction == "downstream" else "Upstream",
    )
    combined_table.add_row("Max Mismatches", f"{args.mismatches}")
    combined_table.add_row("Spacer Lengths", "24")  # Simplified since we're streaming

    # Genome Statistics Section
    combined_table.add_section()
    combined_table.add_row(
        "[bold bright_blue]Genome Statistics[/bold bright_blue]",
        "[italic]Genomic Overview[/italic]",
    )

    # Organism Details
    unique_organisms = set(organisms.values())
    organism_entry = (
        next(iter(unique_organisms))
        if len(unique_organisms) == 1
        else ", ".join(unique_organisms) if unique_organisms else "Unknown"
    )
    combined_table.add_row("Organism", f'"{organism_entry}"')

    for chrom in sorted(seq_lens):
        topology = topologies.get(chrom, "Unknown")

        combined_table.add_row(
            f'Chromosome "{chrom}"', f"{seq_lens[chrom]:,} nt ({topology})"
        )

    # Ambiguity Statistics
    ambiguous_coordinates = {
        (chrom, pos % seq_lens[chrom])
        for chrom, pos in locus_map
        if len(locus_map[(chrom, pos)]) > 1
    }
    ambiguous_locus_tags = {
        entry[0]
        for chrom, pos in ambiguous_coordinates
        for entry in locus_map[(chrom, pos)]
    }

    total_genes = sum(all_genes.values())
    overlapping_genes_percent = (
        (len(ambiguous_locus_tags) / total_genes) * 100 if total_genes > 0 else 0
    )

    combined_table.add_row("Total Genes", f"{total_genes:,}")
    combined_table.add_row(
        "Overlapping Genes",
        f"{len(ambiguous_locus_tags):,} ({overlapping_genes_percent:.2f}%)",
    )
    combined_table.add_row("Ambiguous Coordinates", f"{len(ambiguous_coordinates):,}")

    # Barcode Mapping Statistics
    combined_table.add_section()
    combined_table.add_row(
        "[bold bright_green]Barcode Mapping Statistics[/bold bright_green]",
        "[italic]Targeting Overview[/italic]",
    )

    unique_chrs = results["chr"].nunique()
    unique_tags = results["locus_tag"].nunique()
    total_genes_targeted_percent = (
        (unique_tags / total_genes) * 100 if total_genes > 0 else 0
    )

    combined_table.add_row("Chromosomes Targeted", "Multiple")
    combined_table.add_row("Genes Targeted", "Streaming mode - stats simplified")

    # Detailed Barcode Mapping
    combined_table.add_row("Unique Barcodes", "Streaming mode - stats simplified")
    combined_table.add_row("Off-targeting Barcodes", f"{off_target_spacers:,}")
    combined_table.add_row(
        "Non-targeting Barcodes", "Streaming mode - stats simplified"
    )

    console.log(combined_table)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Map barcodes to a circular genome",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument("sgrna_file", help="Path to sgRNA_fasta_file", type=str)
    parser.add_argument("genome_file", help="Path to genome_gb_file", type=str)
    parser.add_argument("pam", help="PAM sequence", type=str)
    parser.add_argument("mismatches", help="Number of allowed mismatches", type=int)

    parser.add_argument(
        "--mismatches",
        type=int,
        default=1,
        metavar="(0-2)",
        help="Number of mismatches to constitute an offtarget.",
    )
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
