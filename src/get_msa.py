"""
Exemplary run: `python ./src/get_msa.py ./data/clusters/full/ ./data/clusters/full_msa/`
"""
import click
import glob
import multiprocessing
import os

from functools import partial


@click.command()
@click.argument("input_dir")
@click.argument("output_path")
@click.option("-q", "--quiet_mafft", "is_quiet", is_flag=True, required=False,
              help="Set to use `--quiet` option in MAFFT - the program does not report progress.")
@click.option("-m", "--mafft_path", "mafft", default="mafft", help="Command/path to run MAFFT.")
def get_msa(input_dir: str, output_path: str, is_quiet: bool,  mafft: str="mafft"):
    """
    This script runs MAFFT on .fasta files from a given directory (`input_dir`) and saves obtained MSAs
    in `output_path` directory.
    """
    input_files = glob.glob(f"{input_dir}/*.fasta")
    no_cores = int(0.75 * multiprocessing.cpu_count())
    with multiprocessing.Pool(no_cores) as pool:
        pool.map(partial(run_mafft, output_path, is_quiet, mafft), input_files)


def run_mafft(output_path: str, is_quiet: bool, mafft: str, file_name: str):
    """
    Runs MAFFT on given .fasta file and saves it in given directory (`output_name`) as `mafft_{family_name}.fasta`.
    """
    output_name = f"./{output_path}/mafft_{file_name.split('/')[-1]}"
    quiet_param = "--quiet " if is_quiet else ""
    if os.path.dirname(output_name):
        os.makedirs(os.path.dirname(output_name), exist_ok=True)
    command = '{} --auto --inputorder --preservecase {}"{}" > "{}"'.format(mafft, quiet_param, file_name, output_name)
    os.system(command)
    if os.stat(output_name).st_size == 0:
        os.remove(output_name)
        command = '{} --auto --inputorder --preservecase --anysymbol {}"{}" > "{}"'.format(mafft, quiet_param, file_name, output_name)
        os.system(command)
    return output_name


if __name__ == "__main__":
    get_msa()
