#!/bin/bash

# Exemplary run: `./run.sh -s "./data/organisms.txt" -o "./data/" -f bin/fasturec`

helpFunction() {
  echo "Usage: $(basename "$0") [-p PATH] [-o PATH]"
  echo "-s define a path to a file with species names"
  echo "-o define a path output_files"
  echo "-f path to Fasturec executable."
  echo ""
  echo "All requirements can be found in readme.md file"
  exit
}

while getopts "s:o:f:" opt
do
   case "$opt" in
      s ) species_file="$OPTARG" ;;
      o ) output_dir="$OPTARG" ;;
      f ) fasturec_path="$OPTARG" ;;
      ? ) helpFunction ;;
   esac
done

if [ -z "$species_file" ] || [ -z "$output_dir" ] || [ -z "$fasturec_path" ]
then
   echo "Not all required parameters are filled!";
   helpFunction
fi

echo "[0/] Installing Python requirements..."
pip install -r requirements.txt

echo "Downloading proteomes..."
time python ./src/get_proteomes.py "$species_file" "$output_dir"

echo "Creating clusters..."
python ./src/get_clusters.py "$output_dir"/proteomes.fasta "$output_dir"
python ./src/get_filtered_clusters.py "$species_file" "$output_dir"/clusters/parsed/ "$output_dir"/clusters/full/
python ./src/get_filtered_clusters.py "$species_file" "$output_dir"/clusters/parsed/ "$output_dir"/clusters/para/ -t 3

echo "Computing MSA..."
python ./src/get_msa.py "$output_dir"/clusters/full/ "$output_dir"/clusters/full_msa/ -q
python ./src/get_msa.py "$output_dir"/clusters/para/ "$output_dir"/clusters/para_msa/ -q

echo "Creating NJ trees..."
python ./src/get_nj_trees.py "$output_dir"/clusters/full_msa/ "$output_dir"/one2one/
python ./src/get_nj_trees.py "$output_dir"/clusters/para_msa/ "$output_dir"/paralogs/
echo "NJ trees without bootstrap created."
python ./src/get_nj_trees.py "$output_dir"/clusters/full_msa/ "$output_dir"/bootstrap/ -b 70
echo "NJ trees with bootstrap created."
wait

getConsensusTree121() {
  echo "Building consensus tree..."
  Rscript ./src/get_consensus_tree.R "$output_dir"/one2one/nj_trees.nwk "$output_dir"/one2one/
  echo "Consensus tree created."
}

getSupertrees() {
  echo "Building supertree..."
  $fasturec_path -G "$output_dir"/one2one/nj_trees_length_less.nwk -Y
  python ./src/get_supertree.py "./" "$output_dir"/one2one/
  echo "Supertree created."

  echo "Building supertree for paralogs..."
  $fasturec_path -G "$output_dir"/paralogs/nj_trees_length_less.nwk -Y
  python ./src/get_supertree.py "./" "$output_dir"/paralogs/
  echo "Supertree created."

  echo "Building supertree tree with bootstrap..."
  $fasturec_path -G "$output_dir"/bootstrap/nj_trees_length_less.nwk -Y
  python ./src/get_supertree.py "./" "$output_dir"/bootstrap/
  echo "Supertree tree with bootstrap created."
}

getConsensusTreeBootstrap() {
  echo "Building consensus tree with bootstrap..."
  Rscript ./src/get_consensus_tree.R "$output_dir"/bootstrap/nj_trees.nwk "$output_dir"/bootstrap/
  echo "Consensus tree with bootstrap created."
}

getConsensusTree121 & getSupertrees & getConsensusTreeBootstrap
wait

python ./src/get_rf.py "$output_dir"/ref_tree.nwk "$output_dir" -t "$output_dir"/one2one/consensus_tree.nwk -t "$output_dir"/one2one/super_tree.nwk -t "$output_dir"/bootstrap/consensus_tree.nwk -t "$output_dir"/bootstrap/super_tree.nwk -t "$output_dir"/paralogs/super_tree.nwk -t "$output_dir"/timetree/organisms_timetree.nwk
python ./src/get_rf.py "$output_dir"/timetree/organisms_timetree.nwk "$output_dir" -t "$output_dir"/one2one/consensus_tree.nwk -t "$output_dir"/one2one/super_tree.nwk -t "$output_dir"/bootstrap/consensus_tree.nwk -t "$output_dir"/bootstrap/super_tree.nwk -t "$output_dir"/paralogs/super_tree.nwk -t "$output_dir"/ref_tree.nwk

