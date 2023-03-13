"""
Exemplary run: `python ./src/get_proteomes.py ./data/organisms.txt ./data/`
"""
import click
import glob
import os
import subprocess
import time
import wget

from http.client import IncompleteRead


@click.command()
@click.argument("species_file")
@click.argument("output_path")
def get_proteomes(species_file: str, output_path: str):
    """
    The script reads a given file with species names (`species_file`) and downloads their proteomes from UniProt,
    merges them to one file and saves it in a desired directory (`output_path`).
    File format: one line - one name.
    """
    proteomes_dir = f"{output_path}/proteomes/"
    download_proteomes(get_ids(get_species(species_file)), proteomes_dir)
    write_proteomes_file(proteomes_dir, f"{output_path}/proteomes.fasta")


def download_proteomes(ids: list[str], output_path: str):
    """
    Downloads proteomes for given ids from UniProtKB. If not possible tries to proteomes from UniParc.
    """
    if os.path.dirname(output_path):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    for id in ids:
        file_name = f"{output_path}/{id}.fasta.gz"
        try:
            db = download_from_UniProtKB(id, file_name)
        except IncompleteRead:
            continue
        except Exception as e:
            print(e)
            os.remove(file_name)
            db = download_from_UniParc(id, file_name)
        finally:
            decompress(file_name)
            time.sleep(10)

        if os.stat(file_name.split(".gz")[0]).st_size == 0:
            subprocess.run(["rm", f"{file_name.split('.gz')[0]}"])
            db = download_from_UniParc(id, file_name)
            decompress(file_name)
        print(f"Proteome for {id} downloaded from {db}.")
        print("---")
    empty_files = [file_name.split(".")[0] for file_name in glob.glob(f"{output_path}/*.fasta")
                   if os.stat(file_name).st_size == 0]
    if empty_files:
        print(f"Could not receive proteome for: {empty_files}")


def download_from_UniProtKB(id: str, file_name: str) -> str:
    wget.download(f"https://rest.uniprot.org/uniprotkb/stream?compressed=true&download=true&format=fasta&query=%28%28proteome%3A{id}%29%29",
                  out=file_name)
    return "UniProtKB"


def download_from_UniParc(id: str, file_name: str) -> str:
    wget.download(f"https://rest.uniprot.org/uniparc/stream?compressed=true&download=true&format=fasta&query=%28%28upid%3A{id}%29%29",
                  out=file_name)
    return "UniParc"


def decompress(file_name: str):
    os.popen(f"gunzip -f {file_name}")


def get_ids(species: list[str]) -> list[str]:
    """
    Returns list of ids for given ids-species list.
    """
    return [elem.split()[0] for elem in species]


def get_species(species_file: str) -> list[str]:
    """
    Returns list of ids-species acquired from a given file.
    File format: one line - one name.
    """
    with open(species_file, "r") as file:
        return [line.strip() for line in file if line.strip()]


def write_proteomes_file(proteomes_dir: str, output_name: str):
    """
    Saves all proteomes from all .fasta files in a given directory to one .fasta file.
    Adds proteome UniProt ID at the end of every header.
    """
    for file_name in glob.glob(f"{proteomes_dir}/*.fasta"):
        with open(file_name, "r") as proteome_file, open(output_name, "a") as output_file:
            for line in proteome_file:
                if line.strip():
                    if line.startswith(">"):
                        output_file.write(f"{line.strip()} {file_name.split('/')[-1].split('.')[0]}\n")
                    else:
                        output_file.write(line)


if __name__ == "__main__":
    get_proteomes()
