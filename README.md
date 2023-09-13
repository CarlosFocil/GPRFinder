# GPRFinder: A tool for guiding the curation of genome-scale metabolic models through comparative genomics
GPRFinder is a homology-based tool for finding gene-reaction association rules of genes missing in Genome-scale Metabolic Models (GEMs). The tool compares the whole genome of your model (target) versus the genome of a model from the BiGG database (template), and extracts bi-directional best hits (BBHs) with a provided PID threshold. It generates a table in CSV format containing the ID of the orthologous genes (target vs. template) with information on the reactions (Reaction ID, name, subsystem, and string) in which the template genes participate. This table can be used to accelerate gap-filling and curation processes of GEMs.

## Usage
Example on how to run GPRFinder using the _E. coli_ model (iML1515) as a template, with an 80% PID threshold and 4 threads for the BLAST algorithm
```
python GPR_Finder.py \
-template_model iML1515 \
-target_genome genomes/your_genome.gb \
-target_model your_model_file.json \
-pid 80 -num_threads 4
```

## Getting Started

### Installing
Clone the repository on your system:
```
git clone https://github.com/CarlosFocil/GPRtest.git
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
Move your model and its corresponding genome to the GPRFinder directory:
```
mv path/to/your_model_file.json path/to/GPRFinder
mv path/to/your_genome.gb path/to/GPRFinder/genomes

Note: Make sure the genome of your model is located in the genomes directory inside the GPRFinder directory
```

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## References

[1] Thiele, I., & Palsson, B. Ø. (2010). A protocol for generating a high-quality genome-scale metabolic reconstruction. Nature protocols, 5(1), 93-121.

[2] Norsigian, C. J., Fang, X., Seif, Y., Monk, J. M., & Palsson, B. O. (2020). A workflow for generating multi-strain genome-scale metabolic models of prokaryotes. Nature protocols, 15(1), 1-14.
