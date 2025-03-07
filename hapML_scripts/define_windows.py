import pandas as pd
import numpy as np


class GenomicWindow:
    def __init__(self, observed_frequencies, true_frequencies):
        self.observed_frequencies = observed_frequencies
        self.true_frequencies = true_frequencies

    def _get_true_freq_at_position(self, position):
        nearest_pos = min(self.true_frequencies['pos'], key=lambda x: abs(x - position))
        window_true_freq = self.true_frequencies[self.true_frequencies['pos'] == nearest_pos].iloc[:, 2:10]
        return window_true_freq.reset_index(drop=True)

    def define_windows_by_snps(self, snp_number, step, sim):
        results = {}
        start_idx = 0
        end_idx = snp_number

        while end_idx < len(self.observed_frequencies):
            window_data = self.observed_frequencies.iloc[start_idx:end_idx, :]

            # Determine window center
            middle_obs = window_data['pos'].median()
            window_data_trimmed = window_data.iloc[:, 2:11].reset_index(drop=True)

            # Determine true frequencies at window center
            window_true_freq = self._get_true_freq_at_position(middle_obs)

            results[(sim, (window_data['pos'].iloc[0], window_data['pos'].iloc[-1]))] = {
                'window': window_data_trimmed,
                'true_freq_row': window_true_freq
            }

            start_idx += step
            end_idx += step

        return results

    def define_windows_by_size(self, window_size_bp, step_bp, min_snps_per_window, sim):
        results = {}
        pos_min, pos_max = self.observed_frequencies['pos'].min(), self.observed_frequencies['pos'].max()

        for window_start in range(pos_min, pos_max - window_size_bp, step_bp):
            window_end = window_start + window_size_bp
            window_data = self.observed_frequencies[
                (self.observed_frequencies['pos'] >= window_start) & (self.observed_frequencies['pos'] < window_end)
            ]

            if len(window_data) >= min_snps_per_window:
                middle_obs = window_data['pos'].median()
                window_data_trimmed = window_data.iloc[:, 2:11].reset_index(drop=True)
                window_true_freq = self._get_true_freq_at_position(middle_obs)

                results[(sim, (window_start, window_end))] = {
                    'window': window_data_trimmed,
                    'true_freq_row': window_true_freq
                }

        return results


# Example Usage:
# genomic_window = GenomicWindow(observed_frequencies_df, true_frequencies_df)
# results_by_snps = genomic_window.define_windows_by_snps(snp_number=50, step=25, sim="simulation_1")
# results_by_size = genomic_window.define_windows_by_size(window_size_bp=10000, step_bp=5000, min_snps_per_window=20, sim="simulation_1")
