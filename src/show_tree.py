"""
Exemplary run: `python ./data/super_tree.nwk`
"""
import click
import matplotlib.pyplot as plt

from Bio import Phylo


@click.command()
@click.argument("tree_file")
def get_supertree(tree_file: str):
    """
    The script reads and shows tree from .nwk file.
    """
    tree = Phylo.read(tree_file, "newick")
    Phylo.draw(tree, do_show=False)
    plt.show()


if __name__ == "__main__":
    get_supertree()

