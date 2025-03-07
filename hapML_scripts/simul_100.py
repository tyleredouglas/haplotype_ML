import pandas as pd
from genome_simulator import GenomeSimulator
from population_simulator import Simulation


# Load and prepare input data
data = pd.read_csv('/global/scratch/users/tylerdouglas/haplotype_ML/chr3L_RILs_updated.csv')
df_wide = data.pivot_table(index=['CHROM', 'pos'], columns='sample', values='fHap', aggfunc='first')

# Define simulation parameters
n_flies = 300
n_generations = 20
recombination_rate = 0.5
read_num = 500000
num_sims = 100
chromosomes = ['chr3L']
fasta_files = {'chr3L': '/global/scratch/users/tylerdouglas/haplotype_ML/fastas/B.3L.fasta'}

# Initialize Genome Simulator
genome_simulator = GenomeSimulator(fasta_files=fasta_files, chromosomes=chromosomes)

# Run simulations
for i in range(1, num_sims + 1):
    sim_name = f'sim{i}'

    # Simulate population genetics
    population_simulation = Simulation(df_wide, n_flies, n_generations, recombination_rate)
    simulated_population = population_simulation.run()

    # Calculate true haplotype frequencies
    true_freqs = population_simulation.calculate_haplotype_frequencies()

    # Generate paired-end reads and calculate read depth
    read_depth = genome_simulator.generate_paired_end_reads(simulated_population, read_num, sim_name)

    # Save results
    read_depth.to_csv(f'{sim_name}_depth.csv', index=True)
    true_freqs.to_csv(f'{sim_name}_true_freqs.csv', index=True)
