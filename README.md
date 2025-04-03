Overview
This workflow simulates genomic data, generates synthetic paired-end reads, and estimates haplotype frequencies using LightGBM. It covers:

1) Population Simulation: Start with founder genotypes, simulate random mating (with drift and recombination).
2) NGS Read Simulation: Use founder FASTAs to generate paired-end reads, complete with random fragment lengths and errors.
3) Mapping & SNP Calling.
4) Haplotype Frequency Estimation: Define genomic windows, match observed SNP frequencies to “true” haplotype frequencies, then train a LightGBM model to predict those frequencies.

Files
pop_simulator.py: Contains Chromosome, Population, and Simulation classes for simulating populations.
read_simulator.py: Includes GenomeSimulator for simulating reads from the populations.
simul_100.py: Example code to generate training data.
haplotype_ml.ipynb: walkthrough of the above functions in one place.
