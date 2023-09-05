__author__ = "Nadia Baig"
__copyright__ = "Copyright (C) 2023  Nadia Baig"
__license__ = "Public Domain"
__version__ = "1.0"

import pandas as pd
import matplotlib.pyplot as plt

# Load multiple dataframes
dataframes = [pd.read_csv("merged_kmers.txt", sep="\t"),
              pd.read_csv("merged_2_otolia.txt", sep="\t")]
              #pd.read_csv("merged_2_otolia_3.txt", sep="\t")



# Create a list of colors for differentiating dataframes
colors = ['green', 'blue']

# Get unique chromosome names from the first dataframe (assuming they are the same for all dataframes)
chromosomes = dataframes[0]['chromosome'].unique()

# Calculate the number of rows and columns for subplots
num_rows = (len(chromosomes) + 3) // 3  # Ensure at least 3 per row
num_cols = min(len(chromosomes), 3)

# Create subplots
fig, axes = plt.subplots(num_rows, num_cols, figsize=(18, 4 * num_rows), sharex=True)

# Flatten the axes array for easy iteration
axes = axes.flatten()

# Iterate through each chromosome and plot the number of unique Kmers for each dataframe
for i, chromosome in enumerate(chromosomes):
    ax = axes[i]
    for j, dataframe in enumerate(dataframes):
        chromosome_data = dataframe[dataframe['chromosome'] == chromosome]
        ax.plot(chromosome_data['start'], chromosome_data['unique_kmers'], linestyle='-', color=colors[j], label=f'Dataframe {j+1}')
    ax.set_ylabel(f'Number of unique Kmers ({chromosome})')
    ax.set_xlabel('Genome Position')
    ax.legend()

# Remove any empty subplots if the number of chromosomes is not a multiple of 4
for i in range(len(chromosomes), num_rows * num_cols):
    fig.delaxes(axes[i])

# Adjust layout for better spacing
plt.tight_layout()

# Show the plots
plt.show()
