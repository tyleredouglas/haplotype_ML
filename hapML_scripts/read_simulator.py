#simulated paired-end NGS data given template fasta files and a matrix of 'genomes' (dataframe where rows are genomic positions with known genoytpes, 
#each column is an individual genome)

import pandas as pd
import numpy as np 
from Bio import SeqIO
import random

#create dictionary of fasta sequences
def read_fasta(file_path):

    sequences = {}

    for record in SeqIO.parse(file_path, "fasta"):

        #save fastas in dict
        haplotype_id = record.id
        sequences[haplotype_id] = record.seq

    return sequences

founders_3L = read_fasta('/global/scratch/users/tylerdouglas/haplotype_ML/fastas/B.3L.fasta')

#dict so that the correct haplotype can be called from a randomly selected chromosome by read_coordinates
chromosome_dict = {
    'chr3L' : founders_3L
}

chroms = ['chr3L']

#define coordinates for simulated PE read
def read_coordinates(population, chroms):

    #get read length from normal distribution
    read_length_f = np.round((np.random.normal(loc = 150, scale = 10))).astype(int)
    read_length_r = np.round((np.random.normal(loc = 150, scale = 10))).astype(int)

    #library fragment length and inner distance
    fragment = np.round((np.random.normal(loc = 500, scale = 50))).astype(int)
    gap = fragment - (read_length_f + read_length_r)

    #select chromosome to read
    chromosome = np.random.choice(chroms, replace = True) 
    
     #subset population data to get max position of chromosome
    chrom_of_interest = population.loc[chromosome]
    max_pos = chrom_of_interest.index.max()

    #define read boundaries
    read_start_f = np.random.randint(1, max_pos - 2000) 
    read_end_f = read_start_f + read_length_f

    read_end_r = read_end_f + gap
    read_start_r = read_end_r + read_length_r

    return chromosome, read_start_f, read_end_f, read_length_f, read_start_r, read_end_r, read_length_r


#generated PE reads
def generate_reads(population, chromosome_list, read_num, out_name):

    #store number of reads generated per haplotype per position (read depth) 
    haplotype_counts = pd.DataFrame(0, index=population.index, columns = ['B1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7', 'B8'])

    read_count = 0

    #initialize empty fastq file
    with open(f'{out_name}_1.fastq', 'w') as fastq_file1, open(f'{out_name}_2.fastq', 'w') as fastq_file2:

        while read_count < read_num:

            #randomly define read coordaintes
            chrom, start_f, end_f, length_f, start_r, end_r, length_r = read_coordinates(population, chromosome_list)

            #subset population dataframe to chromosome where read is
            chrom_of_interest = population.loc[chrom]

            #pick random chromosome to sequence
            template = chrom_of_interest.iloc[:, random.randint(1, len(chrom_of_interest.columns) - 1)]

            #get haplotype at position nearest to read location
            nearest_pos = min(template.index, key=lambda x: abs(x - start_f))
            haplotype = template[nearest_pos]

            #generate read
            fasta_seqs = chromosome_dict[chrom]
            haplotype_sequence = fasta_seqs[haplotype] 
            
            forward_read = haplotype_sequence[start_f:end_f]
            reverse_read = haplotype_sequence[end_r:start_r] # end_r is start of read before flipping it
            reverse_read = reverse_complement(reverse_read)
            reverse_read = reverse_read[::-1]  # flip reverse the sequence
            

            #track reads per haplotype
            haplotype_counts.at[(chrom, nearest_pos), haplotype] += 1
            read_count += 1
            read_ID = random.randint(1000, 9999) 

            #write read to fastq file:
            read_quality_f = 'I'*len(forward_read)
            fastq_file1.write(f'@{chrom}_{start_f}_{read_ID}/1\n')
            fastq_file1.write(f'{forward_read}\n')
            fastq_file1.write("+\n")
            fastq_file1.write(f'{read_quality_f}\n')

            read_quality_r = 'I'*len(reverse_read) 
            fastq_file2.write(f'@{chrom}_{start_f}_{read_ID}/2\n')  
            fastq_file2.write(f'{reverse_read}\n')
            fastq_file2.write("+\n")
            fastq_file2.write(f'{read_quality_r}\n')
          

    return haplotype_counts, fastq_file1, fastq_file2


#generate reverse complement for PE2
def reverse_complement(seq):

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    return ''.join(complement.get(base, base) for base in seq)


