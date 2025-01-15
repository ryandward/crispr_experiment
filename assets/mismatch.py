import pandas as pd
import numpy as np
import argparse
import logging
import sys
from rich.console import Console
from rich.progress import Progress


def gc_content(seq):
    """Calculate GC content of a DNA sequence."""
    return (seq.count("G") + seq.count("C")) / len(seq)


def calculate_y_pred(original, variant, gc_weight, params):
    """Calculate y_pred based on original and variant sequences."""

    # print(original, variant)

    if original is None or variant is None or original is np.nan or variant is np.nan:
        return None
    elif original == variant:
        return None
    elif len(original) != len(variant):
        return None

    y_pred = params["intercept"]

    for pos, (orig_base, var_base) in enumerate(zip(original, variant)):
        if orig_base != var_base:
            y_pred += params[f"{pos}"]
            y_pred += params[f"{orig_base}{var_base}"]

    y_pred += gc_weight * gc_content(original)
    return y_pred


def read_parameters(file_path):
    """Read parameters from a CSV file."""
    try:
        df = pd.read_csv(file_path)
        params = df.set_index("feature")["weight"].to_dict()
        return params
    except Exception as e:
        logging.error(f"Error reading parameters file: {e}")
        sys.exit(1)


def generate_header():
    """Generate and print header row for output."""
    header = ["original", "variant", "change_description", "y_pred"]
    print("\t".join(header))


def find_closest_mismatch(score, mismatches, mismatch_list):
    """Find the closest mismatch based on the desired score."""
    closest_score = None
    closest_mismatch = None
    for mismatch, mismatch_score in mismatches:
        if closest_score is None or abs(mismatch_score - score) < abs(
            closest_score - score
        ):
            if mismatch not in [x[0] for x in mismatch_list]:
                closest_score = mismatch_score
                closest_mismatch = mismatch
    return closest_mismatch, closest_score


def print_mismatches(mismatch_list, spacer):
    """Print the mismatches in a formatted manner."""
    for i, (mismatch, score) in enumerate(mismatch_list):
        if mismatch is not None:
            target_mismatch = (
                spacer[: mismatch[0]] + mismatch[1] + spacer[mismatch[0] + 1 :]
            )
            change_description = f"{spacer[mismatch[0]]}{mismatch[0]+1}{mismatch[1]}"
            row = [spacer, target_mismatch, change_description, f"{score:.4f}"]
            print("\t".join(row))


def generate_mismatches(
    spacers, min_score, max_score, step, parameters, spacer_original
):
    """Generate mismatches for a list of spacers based on desired scores."""
    nucleotides = ["A", "C", "G", "T"]

    for spacer in spacers:
        desired_scores = np.arange(min_score, max_score + step, step)

        mismatches = []
        for pos in range(len(spacer)):
            for nt in nucleotides:
                if nt == spacer[pos]:
                    continue
                target_mismatch = spacer[:pos] + nt + spacer[pos + 1 :]
                y_pred = calculate_y_pred(
                    spacer, target_mismatch, parameters["GC_content"], parameters
                )
                mismatches.append(((pos, nt), y_pred))

        mismatch_list = []
        for score in desired_scores:
            closest_mismatch, closest_score = find_closest_mismatch(
                score, mismatches, mismatch_list
            )
            if closest_mismatch is not None:
                mismatch_list.append((closest_mismatch, closest_score))

    print_mismatches(
        mismatch_list, spacer_original
    )  # Use spacer_original instead of spacer


def process_guides_output(data, parameters, num_variants):
    """Process the output from find_guides and add mismatch predictions."""
    all_rows = []
    perfect_matches = data[data["mismatches"] == 0].copy()
    nucleotides = ["A", "C", "G", "T"]

    # Calculate step size based on number of desired variants
    step = 1.0 / (num_variants - 1) if num_variants > 1 else 1.0
    desired_scores = np.arange(0, 1.01, step)  # Include 1.0

    for _, row in perfect_matches.iterrows():
        # Add original row
        row_dict = row.to_dict()
        row_dict["change_description"] = None
        row_dict["y_pred"] = None
        all_rows.append(row_dict)

        target = row["target"].upper()

        # Generate all possible 1-mismatches and their scores
        possible_variants = []
        for pos in range(len(target)):
            for nt in nucleotides:
                if nt == target[pos]:
                    continue

                variant = target[:pos] + nt + target[pos + 1 :]
                y_pred = calculate_y_pred(
                    target, variant, parameters["GC_content"], parameters
                )

                if y_pred is not None:
                    possible_variants.append(
                        {
                            "variant": variant,
                            "position": pos,
                            "nucleotide": nt,
                            "y_pred": y_pred,
                        }
                    )

        # Sort variants by y_pred
        possible_variants.sort(key=lambda x: x["y_pred"])

        # Select variants closest to desired scores
        selected_variants = []
        for desired_score in desired_scores:
            if not possible_variants:  # Skip if no more variants available
                continue

            # Find closest variant to desired score
            closest_variant = min(
                possible_variants, key=lambda x: abs(x["y_pred"] - desired_score)
            )

            # Add variant to selected list and remove from possible variants
            selected_variants.append(closest_variant)
            possible_variants.remove(closest_variant)

        # Create new rows for selected variants
        for variant_info in selected_variants:
            new_row = row_dict.copy()
            new_row["spacer"] = variant_info["variant"]
            new_row["mismatches"] = 1
            new_row["change_description"] = (
                f"{target[variant_info['position']]}{variant_info['position']+1}{variant_info['nucleotide']}"
            )
            new_row["y_pred"] = f"{variant_info['y_pred']:.4f}"
            all_rows.append(new_row)

    result_df = pd.DataFrame(all_rows)
    return result_df.sort_values(["target", "mismatches", "y_pred"])


def main(args):
    console = Console(file=sys.stderr)
    console.log("[bold red]Initializing mismatch calculator[/bold red]")

    console.log(f"Reading parameters from {args.parameters_file}...")
    params = read_parameters(args.parameters_file)

    # Read input from stdin
    try:
        data = pd.read_csv(sys.stdin, sep="\t")
    except Exception as e:
        console.log(f"[bold red]Error reading input data: {e}[/bold red]")
        sys.exit(1)

    # Process the data
    result_df = process_guides_output(data, params, args.num_variants)

    # Before writing the output, sort by locus_tag, offset, and y_pred
    # Changed ascending=[True, True, False] to ascending=[True, True, True] to sort y_pred in ascending order
    result_df = result_df.sort_values(
        ["locus_tag", "offset", "y_pred"], ascending=[True, True, True]
    )

    # Print to stdout
    result_df.to_csv(sys.stdout, sep="\t", index=False, na_rep="None")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add mismatch predictions to CRISPR guide data"
    )
    parser.add_argument(
        "--parameters_file",
        required=True,
        help="Path to the parameters file (CSV format).",
    )
    parser.add_argument(
        "--verbosity",
        choices=["debug", "info", "warning", "error", "critical"],
        default="info",
        help="Set the logging verbosity level (default: info).",
    )
    parser.add_argument(
        "--num_variants",
        type=int,
        default=10,
        help="Number of variants to generate per guide (default: 20)",
    )
    args = parser.parse_args()

    main(args)
