__author__ = "Software Authors Name"
__copyright__ = "Copyright (C) 2004 Author Name"
__license__ = "Public Domain"
__version__ = "1.0"

import pandas as pd

# Read files into pandas DataFrames
file_a_path = '/1data/Nadia/10X_Genomics_data/10XGenomics_data/batch2/usftp1.novogene.com/H202SC19122466/raw_data/HMK3HDSXX/Otolia/kmer_analysis/kmer.file.kmer'
file_b_path = '/1data/Nadia/10X_Genomics_data/10XGenomics_data/batch2/usftp1.novogene.com/H202SC19122466/raw_data/HMK3HDSXX/Otolia/kmer_analysis/WU18-3664-HR.WU18-3664-NR.subtract.Saphir.intersect.Deodara.subtract.sgenome.subtract.dump'

file_a_df = pd.read_csv(file_a_path, sep='\t', header=None, names=['Column1', 'Column2'])
file_b_df = pd.read_csv(file_b_path, sep='\t', header=None, names=['Column1', 'Column2'])

# Get unique values of column1 from both files
file_a_values = set(file_a_df['Column1'])
file_b_values = set(file_b_df['Column1'])

# Get File A specific values of column1 only
file_a_only = file_a_values - file_b_values

# Get common values of column1
common_values = file_a_values.intersection(file_b_values)

# Get File B only values of column1
file_b_only = file_b_values - file_a_values

# Save individual files
file_a_only_df = file_a_df[file_a_df['Column1'].isin(file_a_only)]
common_values_df = file_a_df[file_a_df['Column1'].isin(common_values)]
file_b_only_df = file_b_df[file_b_df['Column1'].isin(file_b_only)]

file_a_only_df.to_csv('Otolia_only.txt', sep='\t', index=False, header=None)
common_values_df.to_csv('common_otolia_sensephir.txt', sep='\t', index=False, header=None)
file_b_only_df.to_csv('Sensephir_only.txt', sep='\t', index=False, header=None)
