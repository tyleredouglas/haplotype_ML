#set genomic window to calculate haplotype frequencies for either by number of SNPs per window, or window size (bp)

import pandas as pd
import numpy as np

def define_by_SNPs(observed_frequencies, true_frequencies, snp_number, step, sim):
    
    results = {}
    window_start = 0
    window_end = window_start + snp_number

    while window_end < observed_frequencies.index.max():

        # define  window in observed_frequencies
        window_data = observed_frequencies.iloc[window_start:window_end, :]
     
        # determine window center
        middle_obs = window_data['pos'].median()
        window_data = window_data.iloc[:, 2:11]
        window_data = window_data.reset_index(drop = True)
            
        # determine true frequencies at window center
        middle_true = min(true_frequencies['pos'], key=lambda x: abs(x - middle_obs))
        window_true_freq = true_frequencies[true_frequencies['pos'] == middle_true].iloc[:, 2:10]
        window_true_freq = window_true_freq.reset_index(drop = True)

        # store the window and corresponding true_freqs row in the dictionary
        results[str(sim), (str(observed_frequencies.iloc[window_start]['pos']), 
                 str(observed_frequencies.iloc[window_end]['pos']))] = {'window': window_data, 'true_freq_row': window_true_freq}
        
        window_start += step
        window_end += step

    return results

def define_windows(observed_frequencies, true_frequencies, window, step, sim):

    results = {}

    for window_start in range(observed_frequencies['pos'].min(), observed_frequencies['pos'].max() - window, step):
        # define  window in observed_frequencies
        window_end = window_start + window
        window_data = observed_frequencies[(observed_frequencies['pos'] >= window_start) & (observed_frequencies['pos'] < window_end)]

        #minimum number of SNPs per window
        if window_data.shape[0] > 20:

            # determine window center
            middle_obs = window_data['pos'].median()
            window_data = window_data.iloc[:, 2:11]
            window_data = window_data.reset_index(drop = True)

            # determine true frequencies at window center
            middle_true = min(true_frequencies['pos'], key=lambda x: abs(x - middle_obs))
            window_true_freq = true_frequencies[true_frequencies['pos'] == middle_true].iloc[:, 2:10]
            window_true_freq = window_true_freq.reset_index(drop = True)

            # store the window and corresponding true_freqs row in the dictionary
            results[str(sim), (str(window_start),
                 str(window_end))] = {'window': window_data, 'true_freq_row': window_true_freq}


    return results
