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


def pam_matches(pam_pattern, extracted_pam):
    if not extracted_pam:
        return False
    if pam_pattern == "N" * len(pam_pattern) or not pam_pattern:
        return True
    return bool(re.match(pam_pattern.replace("N", "."), extracted_pam))


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
    with tempfile.TemporaryDirectory() as temp_dir:
        genome_index_temp_name = tempfile.NamedTemporaryFile(
            dir=temp_dir, delete=False
        ).name
        index_prefix = os.path.join(temp_dir, genome_index_temp_name)

        with open(os.devnull, "w") as devnull:
            bowtie_build_command = [
                "bowtie-build",
                topological_fasta_file_name,
                index_prefix,
            ]

            bowtie_build_process = subprocess.Popen(
                bowtie_build_command, stdout=devnull, stderr=devnull
            )
            bowtie_build_process.wait()

            with tempfile.NamedTemporaryFile(delete=False) as bowtie_output_temp_file:
                bowtie_command = [
                    "bowtie",
                    "-S",
                    "-k 100",
                    "--best",
                    "--nomaqround",
                    "-p",
                    str(num_threads),
                    "--tryhard",
                    "-v",
                    str(num_mismatches),
                    "-x",
                    index_prefix,
                    sgrna_fastq_file_name,
                ]

                bowtie_process = subprocess.Popen(
                    bowtie_command,
                    stdout=subprocess.PIPE,
                    # stderr=devnull,
                    universal_newlines=True,
                )

        if bowtie_process.stdout is None:
            raise RuntimeError("Bowtie was unable to start. Check your installation.")

        if bowtie_process.stdout is not None:
            with pysam.AlignmentFile(bowtie_process.stdout, "r") as samfile:
                results = parse_sam_output(
                    samfile,
                    locus_map,
                    topological_fasta_file_name,
                    args.genome_file,
                    pam,
                    pam_direction,
                    regulatory,
                )

            bowtie_process.wait()
            os.remove(bowtie_output_temp_file.name)

    if not results:
        raise RuntimeError(
            "No results were returned from the Bowtie process. Check your input files and parameters."
        )
    return results


def filter_offtargets_by_pam(df):
    targeting_spacers = df[df["target"].notna()]["spacer"].unique()
    return df[~((df["target"].isna()) & (df["spacer"].isin(targeting_spacers)))]


def main(args):
    # Echo all received arguments to stderr for debugging
    sys.stderr.write(f"Received arguments: {args}\n")
    sys.stderr.flush()

    console = Console(file=sys.stderr)
    console.log("[bold red]Initializing barcode target seeker[/bold red]")
    num_threads = cpu_count() // 2
    with tempfile.TemporaryDirectory() as working_dir:
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
        results = pd.DataFrame(results).drop_duplicates()

        # Add feature types to results
        results["feature_types"] = results["locus_tag"].map(
            lambda x: ",".join(sorted(feature_types.get(x, []))) if x else None
        )

        try:
            results = filter_offtargets_by_pam(results)
        except KeyError as e:
            console.log(
                f"[bold red]The following critical attribute is missing for every barcode:[/bold red] {e}"
            )
            console.log(
                f"[bold yellow]Are you sure your PAM should be[/bold yellow] '{args.pam}' '{args.pam_direction}' [bold yellow]of the spacer?[/bold yellow]"
            )
            console.log(
                f"[bold yellow]The genome file[/bold yellow] '{args.genome_file}' [bold yellow]has definitions for these chromosomes:[/bold yellow]"
            )
            console.log(json.dumps(organisms, indent=4))
            sys.exit(1)

        def adjust_min_tar(row):
            if row["tar_start"] > row["tar_end"]:
                return row["tar_start"] - seq_lens[row["chr"]]
            return row["tar_start"]

        results["min_tar"] = results.apply(adjust_min_tar, axis=1)
        results = results.sort_values(by=["chr", "min_tar", "spacer"])

        spacers_seen = (
            results[["name", "spacer"]].drop_duplicates().groupby("spacer").size()
        )
        results = results.drop("name", axis=1).drop_duplicates()
        results.loc[results["target"].notnull(), "site"] = (
            results["chr"].astype(str) + "_" + results["coords"].astype(str)
        )

        # Unique targets per spacer
        site_counts = (
            results.groupby(["spacer", "site", "sp_dir"])
            .ngroup()
            .groupby(results["spacer"])
            .nunique()
        )

        # Unique genes per spacer
        gene_counts = results.groupby("spacer")["locus_tag"].nunique(dropna=True)

        # Unique intergenic sites per spacer
        intergenic_counts = (
            results.loc[results["locus_tag"].isna(), "site"]
            .groupby(results["spacer"])
            .nunique()
        )

        spacer_lengths = set(results["len"])
        spacer_len_range = (
            str(next(iter(spacer_lengths)))
            if len(spacer_lengths) == 1
            else ",".join(str(spacer_len) for spacer_len in sorted(spacer_lengths))
        )

        # Replace note creation with direct statistics
        stats = pd.DataFrame(
            {
                "sites": site_counts,
                "genes": gene_counts,
                "intergenic": intergenic_counts,
            }
        )

        stats = stats.fillna(0).astype(int)
        results = results.merge(stats, left_on="spacer", right_index=True, how="left")

        column_order = ["spacer", "locus_tag", "gene", "feature_types", "chr"]
        if not (results["pam"].isnull().all() or results["pam"].nunique() == 1):
            column_order.append("pam")

        # Remove the mismatches conditional - we always want to show mismatches
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

        final_results = results.reindex(columns=column_order)
        integer_cols = ["mismatches", "offset", "overlap", "tar_start", "tar_end"]
        for col in integer_cols:
            if col in final_results.columns:
                final_results[col] = final_results[col].astype("Int64")

        console.log("Formatting as TSV...")
        final_results.to_csv(sys.stdout, sep="\t", index=False, na_rep="None")

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
    combined_table.add_row("Spacer Lengths", f"{spacer_len_range}")

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

    combined_table.add_row("Chromosomes Targeted", f"{unique_chrs:,}")
    combined_table.add_row(
        "Genes Targeted", f"{unique_tags:,} ({total_genes_targeted_percent:.2f}%)"
    )

    # Detailed Barcode Mapping
    unique_barcodes = results["spacer"].nunique()
    off_target_spacers = (
        results[results["target"].notnull()]
        .groupby("spacer")["coords"]
        .apply(set)
        .apply(len)
        .gt(1)
        .sum()
    )
    nontargeting_spacers = results[results["target"].isnull()]["spacer"].nunique()

    # Mapping Efficiency
    total_barcodes = len(set(results["spacer"]))
    mapping_efficiency = (
        (unique_barcodes / total_barcodes) * 100 if total_barcodes > 0 else 0
    )

    combined_table.add_row("Unique Barcodes", f"{unique_barcodes:,}")
    combined_table.add_row("Off-targeting Barcodes", f"{off_target_spacers:,}")
    combined_table.add_row("Non-targeting Barcodes", f"{nontargeting_spacers:,}")

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
