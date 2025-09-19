from numpy import full
from targets import create_topological_fasta
import argparse
import sys
from rich.console import Console
import os
import tempfile
from multiprocessing import cpu_count
from Bio import SeqIO
import re
import subprocess
import pandas as pd
from io import StringIO


def is_dna(sequence):
    return all(base in "GATC" for base in sequence)


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
        "N": "[ATGC]",  # any nucleotide (using [ATGC] to match existing pattern)
    }

    regex_pattern = ""
    for char in pam_pattern.upper():
        regex_pattern += degenerate_map.get(char, char)

    return regex_pattern


def find_sequences_with_barcode_and_pam(
    topological_fasta_file_name, barcode_length, pam
):
    from rich.console import Console

    console = Console(file=sys.stderr)

    matching_sequences = set()
    pam_regex_pattern = convert_degenerate_to_regex(pam)
    pam_regex = re.compile(pam_regex_pattern)

    console.log(f"Searching for PAM pattern: {pam} -> {pam_regex_pattern}")

    with open(topological_fasta_file_name, "rt") as handle:
        for record_num, record in enumerate(SeqIO.parse(handle, "fasta")):
            console.log(
                f"Processing chromosome {record_num + 1}: {record.id} ({len(record.seq):,} bp)"
            )

            # Consider both the original sequence and its reverse complement
            for strand_name, sequence in [
                ("forward", str(record.seq)),
                ("reverse", str(record.seq.reverse_complement())),
            ]:
                sequence_len = len(sequence)
                positions_to_check = sequence_len - barcode_length - len(pam) + 1

                if strand_name == "forward":
                    console.log(
                        f"Checking {positions_to_check:,} positions on both strands..."
                    )

                for i in range(positions_to_check):
                    # If PAM is downstream
                    if args.pam_direction == "downstream" and pam_regex.match(
                        sequence[i + barcode_length : i + barcode_length + len(pam)]
                    ):
                        spacer = sequence[i : i + barcode_length]
                        if is_dna(spacer):
                            matching_sequences.add(spacer)
                    # If PAM is upstream
                    elif args.pam_direction == "upstream" and pam_regex.match(
                        sequence[i - len(pam) : i]
                    ):
                        spacer = sequence[i : i + barcode_length]
                        if is_dna(spacer):
                            matching_sequences.add(spacer)

    return matching_sequences


# create a sgRNA fasta file such as >sequence\nsequence\n
def create_sgRNA_fasta(matching_sequences, sgRNA_fasta_file_name):
    with open(sgRNA_fasta_file_name, "wt") as handle:
        for sequence in matching_sequences:
            handle.write(f">{sequence}\n{sequence}\n")


def main(args):
    # Echo all received arguments to stderr for debugging
    sys.stderr.write(f"Received arguments: {args}\n")
    sys.stderr.flush()

    console = Console(file=sys.stderr)

    # Use a temp directory with more space than /tmp
    temp_base = (
        os.path.expanduser("~/tmp")
        if os.path.exists(os.path.expanduser("~/tmp"))
        else None
    )
    with tempfile.TemporaryDirectory(dir=temp_base) as working_dir:
        console.log("[bold red]Initializing barcode target builder[/bold red]")

        topological_fasta_file_name = os.path.join(
            working_dir,
            os.path.splitext(os.path.basename(args.genome_file))[0] + ".fasta",
        )

        sgRNA_fasta_file_name = os.path.join(working_dir, "sgRNA.fasta")

        console.log(
            f"[bold blue]Step 1: Creating topological FASTA from {args.genome_file}[/bold blue]"
        )
        create_topological_fasta(args.genome_file, topological_fasta_file_name)
        console.log("[bold green]✓ Topological FASTA created[/bold green]")

        console.log(
            f"[bold blue]Step 2: Finding sequences with PAM '{args.pam}' and length {args.barcode_length}[/bold blue]"
        )
        matching_sequences = find_sequences_with_barcode_and_pam(
            topological_fasta_file_name, args.barcode_length, args.pam
        )
        console.log(
            f"[bold green]✓ Found {len(matching_sequences):,} potential guides[/bold green]"
        )

        console.log("[bold blue]Step 3: Creating sgRNA FASTA file[/bold blue]")
        create_sgRNA_fasta(matching_sequences, sgRNA_fasta_file_name)
        console.log("[bold green]✓ sgRNA FASTA file created[/bold green]")

        console.log(f"Found {len(matching_sequences):,} potential guides in the genome")

        # Safety check for memory issues
        if len(matching_sequences) > 1_000_000:
            console.log(
                f"[bold red]WARNING: {len(matching_sequences):,} guides found![/bold red]"
            )
            console.log(
                "[bold yellow]This may cause memory issues. Consider using:[/bold yellow]"
            )
            console.log("- More specific PAM (e.g., TTTG instead of TTTV)")
            console.log("- Shorter guide length (20 bp instead of 24 bp)")
            console.log("- Fewer mismatches (1-2 instead of 3)")
            console.log("[bold red]Proceeding anyway...[/bold red]")

        console.log(
            f"Stay tuned... running 'targets.py' to find guides for {args.genome_file} with {args.barcode_length} bp barcodes and {args.pam} PAM sequence"
        )

        result = subprocess.run(
            [
                "python",
                "assets/targets.py",
                sgRNA_fasta_file_name,
                args.genome_file,
                args.pam,
                str(args.mismatches),  # Moved mismatches before optional args
                "--pam-direction",
                args.pam_direction,
                "--regulatory",
                str(args.regulatory),
            ],
            stdout=subprocess.PIPE,
            text=True,
            check=True,
        )

        targets = pd.read_csv(StringIO(result.stdout), sep="\t")

        # Removed all filtering code

        # Convert all numeric columns to integers
        targets = targets.apply(
            lambda col: (
                pd.to_numeric(col, errors="coerce").fillna(0).astype(int)
                if col.dtypes != object
                else col
            )
        )

        targets = targets.sort_values(
            ["chr", "tar_start", "tar_end", "locus_tag", "offset", "overlap"]
        )

        # print csv to output
        targets.to_csv(sys.stdout, sep="\t", index=False, na_rep="None")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find guide sequences with specified PAM patterns in genome",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--genome-file", help="Path to genome GenBank file", type=str, required=True
    )
    parser.add_argument(
        "--pam",
        help="PAM sequence (supports IUPAC degenerate bases)",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--barcode-length", help="Length of the guide sequence", type=int, required=True
    )
    parser.add_argument(
        "--mismatches",
        type=int,
        default=1,
        metavar="(0-3)",
        help="Number of mismatches to constitute an off-target.",
    )
    parser.add_argument(
        "--pam-direction",
        choices=["upstream", "downstream"],
        default="downstream",
        help="Direction of the PAM sequence relative to guide",
    )
    parser.add_argument(
        "--regulatory",
        type=int,
        default=0,
        help="Regulatory region size (bp) around target position for gene association.",
    )
    args = parser.parse_args()

    main(args)
