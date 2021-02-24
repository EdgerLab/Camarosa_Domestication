import matplotlib.pyplot as plt
import os
import numpy as np


def graph_barplot_density_differences(
    values, te_type, output_dir, num_bins=20, display=False
):
    n, bins, patches = plt.hist(
        values, num_bins, facecolor="blue", ec="black", alpha=0.5
    )

    # TODO check out broken axis documentation so I can see more of the entries
    # under 100 genes
    plt.ylabel("Number of Genes")
    plt.xlabel("Difference in TE Density Values")
    plt.title(("TE: " + te_type + ", " + "H4 vs Camarosa"))
    plt.savefig(os.path.join(output_dir, (te_type + "_DensityDifferences.png")))
    if display:
        plt.show()
    plt.close()
