import pandas as pd

# Read and process File A
file_a_path = '/1data/Nadia/10X_Genomics_data/10XGenomics_data/batch2/usftp1.novogene.com/H202SC19122466/raw_data/HMK3HDSXX/Otolia/kmer_analysis/kmer.file.kmer'
file_b_path = '/1data/Nadia/10X_Genomics_data/10XGenomics_data/batch2/usftp1.novogene.com/H202SC19122466/raw_data/HMK3HDSXX/Otolia/kmer_analysis/WU18-3664-HR.WU18-3664-NR.subtract.Saphir.intersect.Deodara.subtract.sgenome.subtract.dump'

def process_file(file_path):
    values = set()
    chunk_size = 10000  # Adjust chunk size as needed
    
    for chunk in pd.read_csv(file_path, sep='\t', header=None, names=['Column1', 'Column2'], chunksize=chunk_size):
        values.update(chunk['Column1'])

    return values

file_a_values = process_file(file_a_path)
file_b_values = process_file(file_b_path)

# Get File A specific values of column1 only
file_a_only = file_a_values - file_b_values

# Get common values of column1
common_values = file_a_values.intersection(file_b_values)

# Get File B only values of column1
file_b_only = file_b_values - file_a_values

# Save individual files
def save_individual_file(data, file_path):
    with open(file_path, 'w') as f:
        for value in data:
            f.write(value + '\n')

save_individual_file(file_a_only, 'Otolia_only.txt')
save_individual_file(common_values, 'common_otolia_sen3.txt')
save_individual_file(file_b_only, 'Sen3Sephir_only.txt')


