from numpy import full
from targets_spy import create_topological_fasta
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
from Bio.Seq import Seq


def is_dna(sequence):
    return all(base in "GATC" for base in sequence)


def find_sequences_with_barcode_and_pam(
    topological_fasta_file_name, barcode_length, pam
):
    matching_sequences = set()
    pam_regex = re.compile(pam.replace("N", "[ATGC]"))

    with open(topological_fasta_file_name, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            # Consider both the original sequence and its reverse complement
            for sequence in [str(record.seq), str(record.seq.reverse_complement())]:
                for i in range(len(sequence) - barcode_length - len(pam) + 1):
                    # If PAM is downstream
                    if args.pam_direction == "downstream":
                        potential_pam = sequence[
                            i + barcode_length : i + barcode_length + len(pam)
                        ]
                        if pam_regex.match(potential_pam):
                            full_sequence = sequence[i : i + barcode_length + len(pam)]
                            if is_dna(full_sequence):
                                matching_sequences.add(full_sequence)
                    # If PAM is upstream
                    else:
                        potential_pam = sequence[i - len(pam) : i]
                        if pam_regex.match(potential_pam):
                            full_sequence = sequence[i - len(pam) : i + barcode_length]
                            if is_dna(full_sequence):
                                matching_sequences.add(full_sequence)

    return matching_sequences


# create a sgRNA fasta file such as >sequence\nsequence\n
def create_sgRNA_fasta(matching_sequences, sgRNA_fasta_file_name):
    with open(sgRNA_fasta_file_name, "wt") as handle:
        for sequence in matching_sequences:
            if args.pam_direction == "downstream":
                sequence = str(Seq(sequence).reverse_complement())
            handle.write(f">{sequence}\n{sequence}\n")


def main(args):
    console = Console(file=sys.stderr)

    with tempfile.TemporaryDirectory() as working_dir:
        console.log("[bold red]Initializing barcode target builder[/bold red]")

        topological_fasta_file_name = os.path.join(
            working_dir,
            os.path.splitext(os.path.basename(args.genome_file))[0] + ".fasta",
        )

        sgRNA_fasta_file_name = os.path.join(working_dir, "sgRNA.fasta")

        create_topological_fasta(args.genome_file, topological_fasta_file_name)

        matching_sequences = find_sequences_with_barcode_and_pam(
            topological_fasta_file_name, args.barcode_length, args.pam
        )

        create_sgRNA_fasta(matching_sequences, sgRNA_fasta_file_name)

        console.log(f"Found {len(matching_sequences):,} potential guides in the genome")
        # console.log(
        #     f"Using optimized parameters to find guides for {args.genome_file} with {args.barcode_length} bp barcodes and {args.pam} PAM sequence"
        # )

        cmd = [
            "python",
            "./assets/targets_spy.py",
            sgRNA_fasta_file_name,
            args.genome_file,
            args.pam,
            str(args.mismatches),
            "--pam-direction",
            args.pam_direction,
            "--regulatory",
            str(args.regulatory),
        ]
        console.log(f"Executing command: {' '.join(cmd)}")

        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            text=True,
            check=True,
        )

        targets = pd.read_csv(StringIO(result.stdout), sep="\t")

        if args.pam_direction == "downstream":
            # RevComp results back but maintain PAM trimming
            targets["spacer"] = targets["spacer"].apply(
                lambda x: str(Seq(x).reverse_complement()) if pd.notnull(x) else x
            )
            targets["target"] = targets["target"].apply(
                lambda x: str(Seq(x).reverse_complement()) if pd.notnull(x) else x
            )
            if "pam" in targets.columns:
                targets["pam"] = targets["pam"].apply(
                    lambda x: str(Seq(x).reverse_complement()) if pd.notnull(x) else x
                )
            # Flip F/R in sp_dir
            targets["sp_dir"] = targets["sp_dir"].map({"F": "R", "R": "F"})

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
        description="Map barcodes to a circular genome",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--genome-file", help="Path to genome_gb_file", type=str, required=True
    )
    parser.add_argument("--pam", help="PAM sequence", type=str, required=True)
    parser.add_argument(
        "--barcode-length", help="Length of the barcode", type=int, required=True
    )
    parser.add_argument(
        "--mismatches",
        type=int,
        default=3,
        metavar="(0-3)",
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
