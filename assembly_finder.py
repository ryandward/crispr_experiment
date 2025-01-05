# assembly_finder.py

import datetime
import sys
import gzip
import zlib
import argparse
import configparser
import os
import re
from urllib.request import urlretrieve
from urllib.error import HTTPError, URLError
from Bio import Entrez
from Bio import SeqIO
from rich.progress import Progress, BarColumn, TextColumn
from rich.console import Console
from rich.table import Table

DESCRIPTION_WIDTH = 40
CONFIG_FILE = "config.ini"
DOWNLOAD_LOG_FILE = "download_log.tsv"


def get_email(console):
    """Retrieve or prompt for the user's email for NCBI Entrez API."""
    config = configparser.ConfigParser()
    config.read(CONFIG_FILE)

    if not config.has_section("Entrez"):
        config.add_section("Entrez")

    email = config.get("Entrez", "email", fallback=None)
    if email is None:
        console.print("NCBI requires an email address for its API usage.", style="cyan")
        email = input("Please enter your email address: ")
        config.set("Entrez", "email", email)
        with open(CONFIG_FILE, "w") as configfile:
            config.write(configfile)
    else:
        console.print(f"Using stored email address: {email}", style="cyan")

    return email


def fetch_assembly_data(accession_numbers, source, progress, console, args):
    """Fetch assembly data for the given accession numbers and source."""
    assembly_data = []
    if args.gui:
        print("Fetching assembly accessions...", file=sys.stderr)
    fetch_task = progress.add_task(
        "[1/3] Fetching assembly accessions...", total=len(accession_numbers)
    )

    for accession_number in accession_numbers:
        db = "nuccore" if len(accession_number) == 11 else "assembly"
        try:
            handle = Entrez.esearch(db=db, term=accession_number)
            record = Entrez.read(handle)
            handle.close()
        except (HTTPError, URLError, RuntimeError) as e:
            console.print(f"[bold red]Error searching {accession_number}: {str(e)}")
            progress.update(fetch_task, advance=1, refresh=True)
            continue

        if not record["IdList"]:
            console.print(f"[bold red]Accession number not found: {accession_number}")
            progress.update(fetch_task, advance=1, refresh=True)
            continue

        assembly_id = record["IdList"][0]
        try:
            handle = Entrez.esummary(db="assembly", id=assembly_id)
            summary = Entrez.read(handle)
            handle.close()
        except (HTTPError, URLError, RuntimeError) as e:
            console.print(f"[bold red]Error while processing {assembly_id}: {str(e)}")
            progress.update(fetch_task, advance=1, refresh=True)
            continue

        accession_key = "Genbank" if source == "GenBank" else "RefSeq"
        ftp_path_key = f"FtpPath_{source}"
        synonym_key = "Synonym"

        # Safely extract nested keys
        try:
            synonym = summary["DocumentSummarySet"]["DocumentSummary"][0][synonym_key][
                accession_key
            ]
            ftp_path = summary["DocumentSummarySet"]["DocumentSummary"][0][ftp_path_key]
        except KeyError:
            console.print(f"[bold red]Missing expected data for {assembly_id}.")
            progress.update(fetch_task, advance=1, refresh=True)
            continue

        assembly_data.append(
            {
                "accession": synonym,
                "ftp_path": ftp_path,
            }
        )
        progress.update(fetch_task, advance=1, refresh=True)

    return assembly_data


def download_assemblies_sequences(assembly_data, source, progress, console, args):
    """Download sequences for the given assembly data."""
    if args.gui:
        print("Downloading sequences...", file=sys.stderr)
    assemblies_sequences = {}
    download_task = progress.add_task(
        "[2/3] Downloading sequences...", total=len(assembly_data)
    )

    for assembly in assembly_data:
        ftp_path = assembly["ftp_path"]
        if not ftp_path:
            console.print(
                f"[bold yellow]No FTP path available for {assembly['accession']}. Skipping."
            )
            progress.update(download_task, advance=1, refresh=True)
            continue

        genbank_file_name = ftp_path.split("/")[-1] + f"_genomic.gbff.gz"
        genbank_url = f"{ftp_path}/{genbank_file_name}"

        try:
            compressed_file, _ = urlretrieve(genbank_url)
            with gzip.open(compressed_file, "rt") as compressed_handle:
                sequences = list(SeqIO.parse(compressed_handle, "genbank"))
            assemblies_sequences[assembly["accession"]] = sequences
            os.remove(compressed_file)  # Clean up downloaded compressed file
        except (URLError, zlib.error, ValueError, IOError) as e:
            console.print(
                f"[bold red]Error downloading {assembly['accession']}: {str(e)}"
            )
            progress.update(download_task, advance=1, refresh=True)
            continue

        progress.update(download_task, advance=1, refresh=True)

    return assemblies_sequences


def save_sequences_to_file(assemblies_sequences, progress, console, args):
    """Save the downloaded sequences to individual files."""
    if args.gui:
        print("Saving sequences to files...", file=sys.stderr)
    save_task = progress.add_task(
        "[3/3] Saving sequences to files...", total=len(assemblies_sequences)
    )

    assemblies_info = []

    for accession, sequences in assemblies_sequences.items():
        file_path = os.path.abspath(f"sequences/{accession}.gb")
        try:
            os.makedirs("sequences", exist_ok=True)
            with open(file_path, "w") as output_handle:
                SeqIO.write(sequences, output_handle, "genbank")
            organism = next(
                iter({seq.annotations.get("organism", "Unknown") for seq in sequences}),
                "Unknown",
            )
            assemblies_info.append((accession, organism, file_path))
        except IOError as e:
            console.print(f"[bold red]Error writing file {file_path}: {str(e)}")
            progress.update(save_task, advance=1, refresh=True)
            continue

        progress.update(save_task, advance=1, refresh=True)

    return assemblies_info


def append_to_download_log(assemblies_info, args):
    """Append the downloaded sequences information to the running log file."""
    with open(DOWNLOAD_LOG_FILE, "a") as log_file:
        if os.stat(DOWNLOAD_LOG_FILE).st_size == 0:
            log_file.write("Timestamp\tAccession\tOrganism\tFilePath\tSource\n")

        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        for accession, organism, file_path in assemblies_info:
            source = args.source  # Now using the source from command line args
            log_file.write(
                f"{timestamp}\t{accession}\t{organism}\t{file_path}\t{source}\n"
            )


def display_assemblies_table(assemblies_info, console):
    """Display the contents of the assemblies information in a rich table."""
    table = Table(title="Assemblies Information")
    table.add_column("Accession", justify="left")
    table.add_column("Organism", justify="left")
    table.add_column("File Path", justify="left")

    for accession, organism, file_path in assemblies_info:
        table.add_row(accession, organism, file_path)

    console.print(table)


def main():
    parser = argparse.ArgumentParser(
        description="Download sequences from NCBI by accession numbers."
    )
    parser.add_argument(
        "accession_numbers",
        metavar="ACCESSION",
        nargs="+",
        help="An accession number of a sequence to download.",
    )
    parser.add_argument(
        "-s",
        "--source",
        choices=["GenBank", "RefSeq"],
        default="GenBank",
        help="The source database of the sequences (default: %(default)s).",
    )
    parser.add_argument(
        "--email",
        help="Entrez email address (skips interactive prompt if given)",
        default=None,
    )
    parser.add_argument(
        "--gui",
        action="store_true",
        help="Enable GUI mode (disables rich formatting).",
    )

    args = parser.parse_args()

    # Initialize Console
    if args.gui:
        # Use plain Console without rich formatting
        console = Console(no_color=True, file=sys.stderr)
    else:
        # Use rich formatting
        console = Console(file=sys.stderr)

    if args.email:
        Entrez.email = args.email
    else:
        Entrez.email = get_email(
            console
        )  # Use existing function to read from config or prompt

    try:
        with Progress(
            TextColumn(
                f"[bold cyan]{{task.description:<{DESCRIPTION_WIDTH}}}", justify="left"
            ),
            BarColumn(bar_width=None),
            TextColumn("{task.percentage:>3.0f}%", justify="right"),
            transient=False,
            console=(
                console if not args.gui else Console(no_color=True, file=sys.stderr)
            ),
        ) as progress:
            assembly_data = fetch_assembly_data(
                args.accession_numbers, args.source, progress, console, args
            )
            assemblies_sequences = download_assemblies_sequences(
                assembly_data, args.source, progress, console, args
            )
            assemblies_info = save_sequences_to_file(
                assemblies_sequences, progress, console, args
            )

        append_to_download_log(assemblies_info, args)

        # Print structured data for GUI mode, or use rich table for terminal
        if args.gui:
            # Print header and data rows directly
            print("AssemblyAccession\tOrganism\tFilePath", file=sys.stderr)
            for accession, organism, file_path in assemblies_info:
                print(f"{accession}\t{organism}\t{file_path}", file=sys.stderr)
            print(
                f"Download complete! {len(assemblies_sequences)} sequences added.",
                file=sys.stderr,
            )
        else:
            display_assemblies_table(assemblies_info, console)
            console.print(
                f"[bold green]Download complete! {len(assemblies_sequences)} new sequences added."
            )

    except KeyboardInterrupt:
        console.print("\n[bold red]Interrupted by user. Exiting.")
        sys.exit(1)


if __name__ == "__main__":
    main()
