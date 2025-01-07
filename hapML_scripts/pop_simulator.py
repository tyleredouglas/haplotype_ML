#given a matrix of RILs with known genotypes at intervals, simulates random mating with drift/recombination

import pandas as pd
import numpy as np
import random

def recombine(chrom1, chrom2):
    
    # randomly select a recombination point (window) between the chromosomes at 10,000bp interval
    recombination_point = random.randrange(0, len(chrom1) - 1)  # Choose a point for crossover
    
    # create the offspring chromosome with recombination
    new_chrom1 = pd.concat([chrom1.iloc[:recombination_point], chrom2.iloc[recombination_point:]])
    new_chrom2 = pd.concat([chrom2.iloc[:recombination_point], chrom1.iloc[recombination_point:]])
    
    return new_chrom1, new_chrom2


def simulate_population(RIL_matrix, n_flies, n_generations, recombination_rate):
    # initiate starting population
    population = RIL_matrix.sample(n=n_flies*2, axis=1, replace=True)
    
    for generation in range(n_generations):
        
        # create an empty dataFrame for next generation
        new_population = pd.DataFrame(index = population.index, columns = population.columns)
        
        # simulate mating and recombination
        i = 0
        while i < len(new_population.columns):

            parent1 = population.iloc[:, random.randint(1, len(population.columns) - 1)]
            parent2 = population.iloc[:, random.randint(1, len(population.columns) - 1)] 
            
            if random.random() < recombination_rate:
           
                offspring1, offspring2 = recombine(parent1, parent2)
              
            else:
                # no recombination, copy parents
                offspring1 = parent1.copy()
                offspring2 = parent2.copy()
            
            # add offspring chromosomes to new generation
            new_population.iloc[:, i] = offspring1
            if i + 1 < n_flies:
                new_population.iloc[:, i + 1] = offspring2

            i += 1

        # drift
        population = new_population.sample(n=n_flies*2, axis=1, replace=True)
    
    return population


#calculates haplotype frequencies in simulated population for model training
def get_true_freqs(simulated_pop):

    haplotype_columns = simulated_pop.columns.difference(['sample', 'CHROM', 'pos'])
    haplotype_counts = simulated_pop[haplotype_columns].apply(lambda x: x.value_counts(), axis=1).fillna(0)
    haplotype_frequencies = haplotype_counts.div(haplotype_counts.sum(axis=1), axis=0)

    return haplotype_frequencies 
