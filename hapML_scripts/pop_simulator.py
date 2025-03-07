import pandas as pd
import numpy as np
import random


class Chromosome:
    def __init__(self, genotype):
        self.genotype = genotype

    @staticmethod
    def recombine(chrom1, chrom2):
        recombination_point = random.randrange(0, len(chrom1.genotype) - 1)
        new_genotype1 = pd.concat([
            chrom1.genotype.iloc[:recombination_point],
            chrom2.genotype.iloc[recombination_point:]
        ])
        new_genotype2 = pd.concat([
            chrom2.genotype.iloc[:recombination_point],
            chrom1.genotype.iloc[recombination_point:]
        ])
        return Chromosome(new_genotype1), Chromosome(new_genotype2)


class Population:
    def __init__(self, ril_matrix, population_size):
        self.ril_matrix = ril_matrix
        self.population_size = population_size
        self.chromosomes = [
            Chromosome(ril_matrix.iloc[:, idx])
            for idx in np.random.choice(ril_matrix.shape[1], population_size * 2, replace=True)
        ]

    def simulate_generation(self, recombination_rate):
        new_chromosomes = []
        for _ in range(self.population_size):
            parent1, parent2 = random.sample(self.chromosomes, 2)

            if random.random() < recombination_rate:
                offspring1, offspring2 = Chromosome.recombine(parent1, parent2)
            else:
                offspring1, offspring2 = Chromosome(parent1.genotype.copy()), Chromosome(parent2.genotype.copy())

            new_chromosomes.extend([offspring1, offspring2])

        # Apply genetic drift by random sampling
        self.chromosomes = random.choices(new_chromosomes, k=self.population_size * 2)

    def simulate_generations(self, n_generations, recombination_rate):
        for generation in range(n_generations):
            self.simulate_generation(recombination_rate)

    def get_population_matrix(self):
        return pd.concat([chrom.genotype for chrom in self.chromosomes], axis=1)


class Simulation:
    def __init__(self, ril_matrix, n_flies, n_generations, recombination_rate):
        self.population = Population(ril_matrix, n_flies)
        self.n_generations = n_generations
        self.recombination_rate = recombination_rate

    def run(self):
        self.population.simulate_generations(self.n_generations, self.recombination_rate)
        return self.population.get_population_matrix()

    def calculate_haplotype_frequencies(self):
        simulated_pop = self.population.get_population_matrix()
        haplotype_columns = simulated_pop.columns.difference(['sample', 'CHROM', 'pos'])
        haplotype_counts = simulated_pop[haplotype_columns].apply(lambda x: x.value_counts(), axis=1).fillna(0)
        haplotype_frequencies = haplotype_counts.div(haplotype_counts.sum(axis=1), axis=0)
        return haplotype_frequencies


# Example usage (replace RIL_matrix with actual DataFrame):
# simulation = Simulation(RIL_matrix, n_flies=100, n_generations=10, recombination_rate=0.1)
# final_population = simulation.run()
# haplotype_freqs = simulation.calculate_haplotype_frequencies()
