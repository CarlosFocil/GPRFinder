# GPRFinder: A tool for guiding the curation of genome-scale metabolic models through comparative genomics
GPRFinder is a homology-based tool for finding gene-reaction association rules of genes missing in Genome-scale Metabolic Models (GEMs). The tool compares the whole genome of your model (target) versus the genome of a model from the BiGG database (template), and extracts bi-directional best hits (BBHs) with a provided PID threshold. It generates a table in CSV format containing the ID of the orthologous genes (target vs. template) with information on the reactions (Reaction ID, name, subsystem, and string) in which the template genes participate. This table can be used to accelerate gap-filling and curation processes of GEMs.

## Usage
Example on how to run GPRFinder using the E. coli model (iML1515) as a template, with an 80% PID threshold and 4 threads for the BLAST algorithm
```
python GPR_Finder.py \
-template_model iML1515 \
-target_genome genomes/your_genome.gb \
-target_model your_model_file.json \
-pid 80 -num_threads 4
```

## Getting Started

These instructions will get you a copy of the project up and running on your local machine.

### Installing
Clone the repository on your system:
```
git clone https://github.com/CarlosFocil/GPRFinder.git
```
To install the Python requirements, change to the GPRFinder directory and run:
```
pip install -r requirments.txt
```
If you don't have the BLAST+ on your system, make sure to install it with the following commands:

```
sudo apt-get update
sudo apt-get install ncbi-blast+
```
