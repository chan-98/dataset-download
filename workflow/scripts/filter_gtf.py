import argparse
from typing import Dict, List
import re


def parse_options():
    """Parses command-line.

    Returns:
        argparse.Namespace: Command-line options.
    """
    parser = argparse.ArgumentParser(prog="yaml_dumper.py")
    parser.add_argument(
        "-i",
        "--input",
        action="store",
        type=str,
        help="GTF to filter.",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--filters",
        action="store",
        type=str,
        help=(
            "Dictionary of filters to apply to the GTF file. The dictionary must be"
            "enclosed in simple quotes"
        ),
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        type=str,
        help="File to output the filter GTF into.",
        required=True,
    )
    args = parser.parse_args()
    return args


def get_filters(filters: str) -> Dict:
    """Gets the requested filters into a dictionary

    Parameters:
        - filters: the filters passed in the parameters

    Returns:
        A dictionary of patterns to look out for
        - "in" patterns of things to keep
        - "out patterns for things to leave out
    """
    comment_pattern = "(^#)"
    flts_in = [comment_pattern]
    flts_out = []
    for flt in filters.split("flt_"):
        if len(flt) > 0:
            filter = flt.split(":")
            flt_name = filter[0]
            for v in filter[1].split(","):
                if v.startswith("+"):
                    flts_in.append(f'({flt_name} "{v[1:]}")')
                elif v.startswith("-"):
                    flts_out.append(f'({flt_name} "{v[1:]}")')
                else:
                    raise ValueError(
                        "Filter values should start with '+' or '-' but filter "
                        f"'{flt_name}' has a value starting with '{v[0]}'"
                    )
    to_return = {"in": flts_in}
    if len(flts_out) > 0:
        to_return["out"] = flts_out
    return to_return


def filter_gtf(fin: str, fout: str, filters: Dict[str, List[str]]):
    """Filter the gtf `fin` by the filters in `fitlers` and output it in `fout`

    Parameters:
        - fin: The GTF to filter
        - fout: The name of the output file
        - filters: A dicionary containing the filters to apply
    """
    flt_in = None
    # If more than one pattern are found in "in" we want to join them by '|'
    if len(filters["in"]) == 1:
        no_positive_flt = True
        flt_in = re.compile(f"{filters['in'][0]}")
    else:
        no_positive_flt = False
        flt_in = re.compile("|".join(filters["in"]))
    flt_out = None
    # Add the negative filters to flt_out if they exist
    if "out" in filters:
        if len(filters["out"]) == 1:
            flt_out = re.compile(f"{filters['out'][0]}")
        elif len(filters["out"]) > 1:
            flt_out = re.compile("|".join(filters["out"]))
    with open(fin, "r") as in_file:
        with open(fout, "w") as out_file:
            for line in in_file.readlines():
                if no_positive_flt:
                    # For each line if there is no positive filters we try to keep all
                    # lines and apply the negative filters.
                    keep_line = True
                else:
                    # If there are positive filters we don't keep the lines by default
                    # and only print them if they are added by a positive filter
                    keep_line = False
                if flt_in.search(line):
                    keep_line = True
                if flt_out is not None and keep_line:
                    if flt_out.search(line):
                        keep_line = False
                if keep_line:
                    out_file.write(line)


def main():
    args = parse_options()
    filters = get_filters(args.filters)
    filter_gtf(fin=args.input, fout=args.output, filters=filters)


if __name__ == "__main__":
    main()
