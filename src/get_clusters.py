"""
Exemplary run: `python ./src/get_clusters.py ./data/proteomes.fasta ./data/`
"""
import click
import os


@click.command()
@click.argument("proteomes_file")
@click.argument("output_path")
def get_clusters(proteomes_file: str, output_path: str):
    """
    The script runs clustering with MMSeqs2 on a given .fasta file (`proteomes_file`) and saves every received cluster
    in a separate .fasta file. All output files are stored in `output_path`.
    """
    run_mmseqs(proteomes_file, output_path)
    clusters = parse_clustering_output(output_path)
    write_clusters_files(clusters, output_path)


def run_mmseqs(proteomes_file: str, output_path: str):
    """
    Runs clustering with MMSeqs2's easy-cluster on a given .fasta file (`proteomes_file`) and saves the output
    in `clusters/` folder in a given directory (`output_path`).
    """
    if os.path.dirname(f"{output_path}/clusters/"):
        os.makedirs(os.path.dirname(f"{output_path}/clusters/"), exist_ok=True)
    os.system(f"mmseqs easy-cluster {proteomes_file} {output_path}/clusters/clusters {output_path}/clusters/clusters_tmp")


def parse_clustering_output(output_path: str) -> dict[str, dict[str, str]]:
    """
    Parses one of the MMSeqs2's output files (the one with `_all_seqs.fasta` suffix, containing all input sequences
    ordered by clusters) to a dictionary ({cluster_name: {header: sequence}}).
    """
    with open(f"{output_path}/clusters/clusters_all_seqs.fasta", "r") as clustering_output:
        clusters = {}
        cluster_name = None
        header = None
        for line in clustering_output:
            if line.startswith(">"):
                if header:
                    cluster_name = header
                    clusters[cluster_name] = {}
                header = line.strip()
            else:
                if header:
                    clusters[cluster_name][header] = line.strip()
                    header = None
                else:
                    raise Exception(f"Wrong file format - sequence has no header defined! Error for line: {line.strip}")
    return clusters


def write_clusters_files(clusters: dict[str, dict[str, str]], output_path: str):
    """
    Writes a .fasta file in `clusters/parsed/` folder in a given directory (`output_path`) for every cluster
    from a given dictionary (`clusters`: {cluster_name: {header: sequence}}).
    """
    output_dir = f"{output_path}/clusters/parsed/"
    if os.path.dirname(output_dir):
        os.makedirs(os.path.dirname(output_dir), exist_ok=True)
    for cluster, headers in clusters.items():
        with open(f"{output_dir}{len(headers)}_{cluster.strip('>')}.fasta", "a") as cluster_file:
            for header, sequence in headers.items():
                cluster_file.write(f"{header}\n{sequence}\n")


if __name__ == "__main__":
    get_clusters()
