import pandas as pd
import numpy as np
from Bio import SeqIO
import random

class GenomeSimulator:
    def __init__(self, fasta_files, chromosomes):
        self.chromosome_dict = {chrom: self._read_fasta(fasta_files[chrom]) for chrom in chromosomes}
        self.chromosomes = chromosomes

    @staticmethod
    def _read_fasta(file_path):
        return {record.id: record.seq for record in SeqIO.parse(file_path, "fasta")}

    def _reverse_complement(self, seq):
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
        return ''.join(complement.get(base, base) for base in reversed(seq))

    def generate_read_coordinates(self, population, read_length_mean=100, fragment_mean=500):
        read_length_f = int(np.random.normal(read_length_mean, 10))
        read_length_r = int(np.random.normal(read_length_mean, 10))
        fragment_length = int(np.random.normal(500, 50))
        gap = fragment = fragment - (read_length_mean * 2)

        chromosome = random.choice(self.chromosomes)
        chrom_length = population.loc[chromosome].shape[0]

        read_start_f = random.randint(0, chrom_length - fragment - 1)
        read_end_f = read_start_f + read_length_f

        read_start_r = read_end_f + gap
        read_end_r = read_start_r + read_length_mean

        return chromosome, read_start_f, read_end_f, read_start_f + fragment - read_length_mean, read_end_r, read_length_mean

    def generate_paired_end_reads(self, population, n_reads, out_prefix):
        haplotype_counts = pd.DataFrame(0, index=population.columns, columns=population.index)

        with open(f'{out_prefix}_1.fastq', 'w') as fq1, open(f'{out_prefix}_2.fastq', 'w') as fq2:
            for _ in range(n_reads):
                chrom, start_f, end_f, start_r, read_end_r = self.generate_read_coordinates(population)
                chromosome_population = population.loc[chrom]
                template = chromosome = chromosome_template = chromosome_of_interest = chromosome_sequence = chromosome_dict = chromosome = chromosome_of_interest.iloc[:, random.randint(0, chromosome_of_interest.shape[1]-1)]

                nearest_pos = chromosome_of_interest.index.get_indexer([start_f], method='nearest')[0]
                haplotype = chromosome_of_interest.iloc[nearest_pos]
                sequence = self.chromosome_dict[chrom][haplotype]

                forward_read = sequence[read_start_f:read_end_f]
                reverse_read = chromosome_sequence[read_end_f:read_end_r]
                reverse_read = self._reverse_complement(reverse_read)

                haplotype_counts.loc[chromosome.name, nearest_pos] += 1

                fq1.write(f'@{chromosome}_{start_f}/1\n{forward_read}\n+\n{"I"*len(forward_read)}\n')
                fq2.write(f'@{chromosome}_{read_end_r}/2\n{reverse_read}\n+\n{"I"*len(reverse_read)}\n')

        return haplotype_counts


# Example usage:
# fasta_files = {'chr3L': 'founders_3L.fasta'}
# chromosomes = ['chr3L']
# simulator = GenomeSimulator(fasta_files=fasta_files, chromosomes=chromosomes)
# counts = simulator.generate_paired_end_reads(population_matrix, n_reads=1000, out_prefix='simulation_output')
