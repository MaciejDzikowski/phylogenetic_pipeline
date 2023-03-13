"""
Exemplary run: `python ./src/get_nj_trees.py ./data/clusters/full_msa/ ./data/`
"""
import Bio.Align
import click
import glob
import multiprocessing
import os

from Bio import Align, AlignIO, Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from Bio.Phylo.Consensus import bootstrap_trees, get_support
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from functools import partial
from typing import Callable


@click.command()
@click.argument("input_dir")
@click.argument("output_dir")
@click.option("-b", "--bootstrap", "bootstrap", required=False, type=click.FloatRange(0, 100),
              help="Set to use bootstrap. The value is a threshold - trees with average support lower than it won't be saved.")
def get_nj_trees(input_dir: str, output_dir: str, bootstrap: float):
    """
    This script takes .fasta files with MSAs from a `input_dir` and saves created trees in `output_dir`.
    """
    input_files = glob.glob(f"{input_dir}/*.fasta")
    if os.path.dirname(f"{output_dir}/nj_trees/"):
        os.makedirs(os.path.dirname(f"{output_dir}/nj_trees/"), exist_ok=True)
    no_cores = int(0.75*multiprocessing.cpu_count())
    with multiprocessing.Pool(no_cores) as pool:
        pool.map(partial(save_tree, get_nj_tree, output_dir, bootstrap), input_files)
    output_name = f"{output_dir}/nj_trees.nwk"
    write_trees_file(f"{output_dir}/nj_trees/", output_name)
    write_length_less_trees_file(output_name, output_dir)


def save_tree(tree_func: Callable, output_name: str, bootstrap: float, mafft_file: str):
    """
    Saves tree created with `tree_func` using an alignment from `mafft_file` as `output_name`.
    Trees with any branch with a negative length value are not saved.
    """
    output_name = f"{output_name}/nj_trees/{mafft_file.split('/')[-1].split('.')[0]}.nwk"
    alignment = read_mafft_output(mafft_file)
    tree = tree_func(alignment, bootstrap)
    if tree and check_nj_tree(tree):
        Phylo.write(tree, output_name, "newick")


def read_mafft_output(mafft_file: str) -> Align.MultipleSeqAlignment:
    """
    Reads an alignment from a given file (`mafft_file`).
    """
    alignment = AlignIO.read(mafft_file, "fasta")
    return alignment


def check_nj_tree(tree: Tree) -> bool:
    """
    Checks if any of the tree's branches has a negative length value. Returns True if not.
    """
    return check_nj_clade(tree.clade)


def check_nj_clade(clade: Clade) -> bool:
    """
    Checks if any of a given clade's branches has a negative length value. Returns True if not.
    """
    clades = clade.clades
    if clade.branch_length < 0:
        return False
    elif clades:
        return sum([check_nj_clade(c) for c in clades]) / len(clades) == 1
    else:
        return True


def get_nj_tree(alignment: Bio.Align.MultipleSeqAlignment, bootstrap: float) -> Tree:
    """
    Returns tree for a given MSA (`alignment`).
    """
    constructor = DistanceTreeConstructor()
    calculator = DistanceCalculator('identity')
    if not bootstrap or check_bootstrap(alignment, bootstrap):
        return constructor.nj(calculator.get_distance(alignment))


def check_bootstrap(alignment: Bio.Align.MultipleSeqAlignment, bootstrap: float) -> bool:
    """
    Generates 100 bootstrap trees for given alignment and returns True if average support value for then is
    higher or equal given `bootstrap` value.
    """
    calculator = DistanceCalculator('identity')
    constructor = DistanceTreeConstructor(calculator)
    trees = list(bootstrap_trees(alignment, 100, constructor))
    avg_supps = [get_avg_supp(get_support(tree, trees)) for tree in trees]
    if sum(avg_supps) / len(trees) >= bootstrap:
        return True


def get_avg_supp(tree: Tree) -> float:
    """
    Returns average support for a given tree.
    """
    supp, no_clades = get_clade_support(tree.clade)
    return supp / no_clades


def get_clade_support(clade: Clade) -> tuple[float, int]:
    """
    Returns sum of support values for clade and its subclades and number of the clades.
    """
    clades = clade.clades
    if clades:
        supp = clade.confidence
        no_clades = 1
        for c in clades:
            s, n = get_clade_support(c)
            supp += s
            no_clades += n
        return supp, no_clades
    else:
        return 0, 0


def write_trees_file(trees_dir: str, output_name: str):
    """
    Saves all trees from all .nwk files in a given directory to one .nwk file.
    """
    tree_files = glob.glob(f"{trees_dir}/*.nwk")
    with open(output_name, "a") as output_file:
        for tree_file in tree_files:
            with open(tree_file, "r") as input_file:
                for line in input_file:
                    output_file.write("".join(line.split("$")))  # remove all `$` used to differentiate headers


def write_length_less_trees_file(trees_file: str, output_dir: str):
    """
    Modifies given .nwk file to a format required for Fasturec - without branch lengths and Inner nodes.
    """
    tmp_name = f"{output_dir}/tmp_nj_trees_length_less.nwk"
    output_name = f"{output_dir}/nj_trees_length_less.nwk"
    os.system(f"perl -ne '$_=~s/:[\d\.]+//g; print $_;' {output_dir}/nj_trees.nwk > {tmp_name}")
    with open(tmp_name, "r") as tmp_file, open(output_name, "w") as output_file:
        for line in tmp_file:
            line_parts = line.strip().strip(";").split("Inner")
            output_file.write(f"{''.join([part.strip('1234567890') for part in line_parts])}\n")
    os.remove(tmp_name)


if __name__ == "__main__":
    get_nj_trees()
