"""
Exemplary run: `python ./src/get_rf.py ./data/timetree/organisms_timetree.nwk ./data/ -t ./data/one2one/consensus_tree.nwk`
"""
import click

from ete3 import Tree


@click.command()
@click.argument("original_tree")
@click.argument("output_dir")
@click.option("-t", "trees", multiple=True, help="Path to a .nwk tree to compare with `original_tree`.")
def get_rf(original_tree: str, output_dir: str, trees: tuple):
    """
    This script computes Robinson-Folds distance between `original_tree` and each .nwk trees from a given `trees` tuple
    and saves the output in `output_dir/report.txt`.
    """
    with open(f"{output_dir}/report_{original_tree.split('/')[-1].split('.')[0]}.txt", "w") as output_file:
        original_obj = Tree(original_tree, format=1)
        for node in original_obj:
            if not node.is_leaf():
                node.support = float(node.name.split("/")[1])
        for tree in trees:
            tree_obj = Tree(tree, format=1)
            for node in tree_obj:
                if not node.is_leaf():
                    node.support = float(node.name.split("/")[1])
            rf_output = original_obj.robinson_foulds(tree_obj, unrooted_trees=True)
            output_file.write(f"{tree} - {original_tree}\n")
            output_file.write(f"RF: {rf_output[0]}\n")
            output_file.write(f"RF max: {rf_output[1]}\n")
            output_file.write("-----\n")


if __name__ == "__main__":
    get_rf()
