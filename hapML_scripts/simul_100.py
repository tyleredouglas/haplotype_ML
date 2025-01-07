import pandas as pd
import numpy as np
import pysam
from Bio import SeqIO
import random
import itertools
from pop_simulator import recombine, simulate_population, get_true_freqs
from read_simulator import read_fasta, read_coordinates, generate_reads, reverse_complement


data = pd.read_csv('/global/scratch/users/tylerdouglas/haplotype_ML/chr3L_RILs_updated.csv')
df = pd.DataFrame(data)
df_wide = df.pivot_table(index=['CHROM', 'pos'], columns='sample', values='fHap', aggfunc='first')

#store founder fastas for each chromosome 
founders_3L = read_fasta('/global/scratch/users/tylerdouglas/haplotype_ML/fastas/B.3L.fasta')

#list of chromosomes to sample (including only 3L for now)
chroms = ['chr3L']

#dict so that the correct haplotype can be called from a randomly selected chromosome by read_coordinates
chromosome_dict = {
    'chr3L' : founders_3L
}

#simulation parameters
n_flies = 300  # population size
n_generations = 20  # number of generations
recombination_rate = 0.5  # probability of recombination occurring
read_num = 500000
num_sims = 100 # number of populations to simulate

#simulate populations, generate reads
for i in range(1, num_sims + 1):

    sim_name = f'sim{i}'
    sim_pop = simulate_population(df_wide, n_flies, n_generations, recombination_rate)
    true_freqs = get_true_freqs(sim_pop)
    read_depth, fq1, fq2 = generate_reads(sim_pop, chroms, read_num, sim_name)
    read_depth.to_csv(f'{sim_name}_depth.csv', index = True)
    true_freqs.to_csv(f'{sim_name}_true_freqs.csv', index = True)


