# Phylogenetic pipeline
Pipeline to build phylogenetic trees for given species whose proteomes can be found in *UniProt* database.

The result of the pipeline:
- consensus tree and supertree with and without bootstrap for NJ trees built on 1-1 clusters (clusters 
with exactly one gene per genome),
- supertree for paralogs (built on all obtained clusters with at least 3 sequences).

Each script included in the program can be used separately. Exemplary usage has been placed at the beginning of each.  

### What do you need?
- `run.sh` script.
- `src/` directory with all provided scripts.
- Text file with organisms from *UniProt* database.

  - Each line should have the following form: `<Entry> <Organism>`.
  - All proteomes should be available in *UniProtKB* or *UniParc* database.


### Requirements
All you need is *Python* (>=3.9.0), *MMseqs2* (recommended: 14-7e284), *MAFFT* (recommended: v7.490), 
*R* (recommended: *4.2.2*), *perl* (recommended: v5.34.0) and [*Fasturec*](https://bio.tools/fasturec).  
Every required library for both languages will be installed automatically. 

For `ape` *R* package it can be helpful to install also:
- `apt-get install liblapack-dev libblas-dev gfortran`


### Exemplary run
- Allow executing `run.sh` script as a program:  
  `chmod +x ./run.sh` 

- Execute the pipeline:  
  `./run.sh -s <organisms_text_file> -o <output_directory> -f <path_to_fasturec>`
