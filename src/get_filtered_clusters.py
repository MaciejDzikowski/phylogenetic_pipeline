"""
Exemplary run: `python ./src/get_filtered_clusters.py ./data/organisms.txt ./data/clusters/parsed/ ./data/clusters/full/ <-t 3>`
"""
import click
import glob
import os


@click.command()
@click.argument("species_file")
@click.argument("input_dir")
@click.argument("output_dir")
@click.option("-t", "not121", required=False, type=int,
              help='Set to allow not only "1 to 1" clusters. The value is a threshold - minimal number of seuences in a cluster.')
def get_filtered_clusters(species_file: str, input_dir: str, output_dir: str, not121: int):
    """
    This script filters out all clusters from `input_dir` with one sequence for each species ("1 to 1") from `species_file`.
    Then, replaces all headers with corresponding species names and saves the output in `output_path` directory.
    """
    species = parse_species(species_file)
    if not121:
        filtered = get_all_clusters(input_dir, species, not121)
    else:
        filtered = get_not121_clusters(input_dir, species)
    unify_headers(output_dir, species, filtered)


def parse_species(species_file: str) -> dict[str, str]:
    """
    Parses a given `species_file` to a dictionary ({uniprot_id: species_name}).
    """
    species = {}
    with open(species_file, "r") as file:
        for line in file:
            line_data = line.strip().split()
            uniprot_id, species_name, strain = line_data[0], "_".join(line_data[1:3]), line_data[3]
            species[uniprot_id] = species_name
    return species


def get_not121_clusters(input_dir: str, species: dict[str, str]) -> list[str]:
    """
    Filters out all clusters from `input_dir` with exactly one sequence for each species from `species` dictionary.
    """
    no_species = len(species)
    clusters_files = glob.glob(f"{input_dir}/{no_species}_*.fasta")
    filtered = []
    for cluster_file in clusters_files:
        is_full = set()
        with open(cluster_file, "r") as file:
            for line in file:
                if line.startswith(">"):
                    is_full.add(species[line.split()[-1]])
        if len(is_full) == no_species:
            filtered.append(cluster_file)
    return filtered


def get_all_clusters(input_dir: str, species: dict[str, str], threshold: int) -> list[str]:
    """
    Filters out all clusters from `input_dir` with exactly one sequence for each species from `species` dictionary.
    """
    no_species = len(species)
    clusters_files = glob.glob(f"{input_dir}/*.fasta")
    for i in range(threshold):
        clusters_files = [file for file in clusters_files if not file.split("/")[-1].startswith(f"{str(i)}_")]
    return clusters_files


def unify_headers(output_dir: str, species: dict[str, str], files: list[str]):
    """
    For each file from a `files` list, replaces all headers with corresponding species names from `species` dictionary
    (using proteome UniProt ID from the end of each header) and saves the output in `output_path` directory.
    Adds "$" at the end of header for every repeating header (if `not121` parameter is set).
    """
    if os.path.dirname(output_dir):
        os.makedirs(os.path.dirname(output_dir), exist_ok=True)
    for file in files:
        with open(file, "r") as input_file, open(f"{output_dir}/unified_{file.split('/')[-1]}", "w") as output_file:
            headers = {}  # {species_name: no_occurrences}
            for line in input_file:
                if line.startswith(">"):
                    species_name = species[line.split()[-1]]
                    if species_name in headers:
                        headers[species_name] += 1
                    else:
                        headers[species_name] = 1
                    distinction_mark = "$" * (headers[species_name] - 1)  # empty string if `not121` parameter not specified
                    output_file.write(f">{species_name}{distinction_mark}\n")
                else:
                    output_file.write(line)
    # FIX should have checked if sequences from different species

if __name__ == "__main__":
    get_filtered_clusters()
