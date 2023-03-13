"""
Exemplary run: `python ./src/get_supertree.py ./ ./data/`
"""
import os

import click
import glob


@click.command()
@click.argument("fasturec_dir")
@click.argument("output_dir")
def get_supertree(fasturec_dir: str, output_dir: str):
    """
    The script takes the most optimal supertree from Fasturec output and saves in `output_dir`.
    """
    fasturec_output = glob.glob(f"{fasturec_dir}/*.fu.txt")
    if fasturec_output:
        fasturec_name = fasturec_output[-1]
        output_name = f"{output_dir}/super_tree.nwk"
        with open(fasturec_name, "r") as fasturec_file, open(output_name, "w") as output_file:
            for line in fasturec_file:
                output_file.write(f"{''.join(line.strip().split()[1:])};")
                break
    else:
        raise Exception("Fasturec output not found!")
    os.remove(fasturec_name)


if __name__ == "__main__":
    get_supertree()
