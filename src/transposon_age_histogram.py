__author__ = "Scott Teresi"

"""
Calculate the (mean) age of a transposon family.
Generate a histogram of the age distribution of a TE family.
User provides the TE family name and the path to the TE annotation file.
"""

import argparse
import os
import logging
import coloredlogs
import pandas as pd
import re

import matplotlib.pyplot as plt


def parse_unfiltered_TE_annotation(TE_annotation_path):
    """
    Parse the unfiltered TE annotation file.
    """
    # Read the TE annotation file
    data = pd.read_csv(
        TE_annotation_path,
        sep="\t",
        comment="#",
        names=[
            "Chromosome",
            "Software",
            "Shorthand_Name",
            "Start",
            "End",
            "Score",
            "Strand",
            "Phase",
            "Attribute",
        ],
        dtype={
            "Chromosome": str,
            "Software": str,
            "Shorthand_Name": str,
            "Start": int,
            "End": int,
            "Score": str,
            "Strand": str,
            "Phase": str,
            "Attribute": str,
        },
    )
    data.drop(columns=["Software", "Score", "Strand", "Phase"], inplace=True)
    # Return the TE annotation
    data["ID"] = data["Attribute"].str.extract(r"ID=(.*?);")
    data["Family"] = data["Attribute"].str.extract(r"Name=(.*?);")
    data["FullType"] = data["Attribute"].str.extract(
        r"classification=(.*?);", flags=re.IGNORECASE
    )
    data["Identity"] = data["Attribute"].str.extract(
        r"identity=(.*?);", flags=re.IGNORECASE
    )
    data.drop(columns=["Attribute"], inplace=True)
    return data


if __name__ == "__main__":
    path_main = os.path.abspath(__file__)
    dir_main = os.path.dirname(path_main)
    parser = argparse.ArgumentParser(description="TODO")

    parser.add_argument(
        "unfiltered_TE_annotation",
        type=str,
        help="Parent path of the cleaned gene annotation file",
    )
    parser.add_argument(
        "family_of_interest",
        type=str,
        help="""string of the TE family, case sensitive, put between quotes""",
    )
    parser.add_argument("output_directory", type=str)

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="set debugging level to DEBUG"
    )

    args = parser.parse_args()
    args.unfiltered_TE_annotation = os.path.abspath(args.unfiltered_TE_annotation)
    args.output_directory = os.path.abspath(args.output_directory)
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logger = logging.getLogger(__name__)
    coloredlogs.install(level=log_level)
    # -------------------------------------------------------------- #

    data = parse_unfiltered_TE_annotation(args.unfiltered_TE_annotation)

    # Filter by shorthand name to be only 'repeat_region' because we don't want
    # to double count anything
    data = data.loc[data["Shorthand_Name"] == "repeat_region"]

    # Print out our family of interest
    subsetted_data = data.loc[data["Family"] == args.family_of_interest]

    # Calculate the mean identity of the family, and log to console
    try:
        subsetted_data["Identity"] = pd.to_numeric(subsetted_data["Identity"])
        mean_identity = subsetted_data["Identity"].mean()
    except ValueError:
        logger.error(
            f"""Could not calculate mean identity of family
            {args.family_of_interest}, likely because the identity field for
            this particular TE contains non-numeric"""
        )
    logger.info(f"Mean identity of family {args.family_of_interest} is {mean_identity}")

    # Calculate the number of elements in the family, so we may add that to the
    # figure
    num_instances = len(subsetted_data)
    # Create artist object
    num_instances_artist = plt.Line2D(
        [], [], linestyle="", label=f"Number of Instances: {num_instances}"
    )

    # Generate a histogram of the identity distribution
    subsetted_data.hist(column="Identity", bins=10, legend=False)
    plt.xlabel(f"Identity Value")
    plt.ylabel(f"Frequency of TE")
    plt.legend(handles=[num_instances_artist], loc="upper left", fontsize="small")

    plt.title(
        f"""Histogram of Identity Values for TE Family: {args.family_of_interest}"""
    )
    file_out = os.path.join(
        args.output_directory, f"{args.family_of_interest}_hist.png"
    )
    logger.info(f"Saving histogram to {file_out}")
    # plt.show()
    plt.savefig(file_out)
