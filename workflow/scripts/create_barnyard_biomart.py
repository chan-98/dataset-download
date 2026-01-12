import argparse
import pandas as pd

"""
This script will take one Human gene annotation BioMart file and one Mouse gene
annotations BioMart to combine them in a single file prefixing the genome 
features with the one expected for the Barnyard genome. See
Human reference (GRCh38) and mouse dataset required for Cell Ranger.
https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38mm10_2020A
"""


def parse_options():
    """Parses command-line.

    Returns:
        argparse.Namespace: Command-line options.
    """
    parser = argparse.ArgumentParser(prog="create_barnyard_biomart.py")
    parser.add_argument(
        "--human",
        action="store",
        type=str,
        help="Human biomart to transform",
        required=True,
    )
    parser.add_argument(
        "--mouse",
        action="store",
        type=str,
        help="Mouse biomart to transform",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output",
        action="store",
        type=str,
        help="File to output the combined biomart files into",
        required=True,
    )
    args = parser.parse_args()
    return args


def generate_biomart(fhuman: str, fmouse: str, fout: str):
    """Generate a new biomart dump from human and mouse biomart annotations

    Parameters:
        - fhuman: Human BioMart annotations
        - fmouse: Mouse BioMart annotations
        - fout: The name of the output file
    """

    biomart_cols = ['ensembl_gene_id', 'ensembl_transcript_id', 'barnyard_gene_name',
                    'ensembl_gene_id_version', 'ensembl_transcript_id_version',
                    'ensembl_peptide_id']

    human_df = pd.read_csv(fhuman, sep='\t')
    human_df = human_df.rename(
        columns={"hgnc_symbol": "barnyard_gene_name", "hgnc_id": "barnyard_gene_id"})
    mouse_df = pd.read_csv(fmouse, sep='\t')
    mouse_df = mouse_df.rename(
        columns={"mgi_symbol": "barnyard_gene_name", "mgi_id": "barnyard_gene_id"})

    # Need to put a prefix for the following columns:
    for column in biomart_cols:
        human_df[column] = 'GRCh38_' + human_df[column]
    for column in biomart_cols:
        mouse_df[column] = 'mm10___' + mouse_df[column]

    # Join the 2 dataframes and save it
    df = pd.concat([human_df, mouse_df])
    df.to_csv(fout, sep='\t', index=False)


def main():
    args = parse_options()
    generate_biomart(fhuman=args.human, fmouse=args.mouse, fout=args.output)


if __name__ == "__main__":
    main()
