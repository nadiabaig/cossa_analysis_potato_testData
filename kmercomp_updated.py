#Code for comparing several genotypes to get unique kmers linked with resistance
__author__ = "Nadia Baig"
__copyright__ = "Copyright (C) 2023  Nadia Baig"
__license__ = "Public Domain"
__version__ = "1.0"

import pandas as pd

#reading all necessary files
# Read and process File A
file_otolia ='Otolia.file.kmer'
file_sen3 = 'WU18-3664-HR.WU18-3664-NR.subtract.Saphir.intersect.Deodara.subtract.sgenome.subtract.dump'
file_kuba='Kuba.file.kmer'
file_laura='Laura.file.kmer'


def process_file(file_path):
    values = set()
    chunk_size = 10000  # Adjust chunk size as needed
    
    for chunk in pd.read_csv(file_path, sep='\t', header=None, names=['Sequence', 'count'], chunksize=chunk_size):
        values.update(chunk['Sequence'])

    return values

file_otolia_values = process_file(file_otolia)
file_sen3_values = process_file(file_sen3)
print(type(file_otolia_values))

# Find common sequences b/w otolia and sen3
common_sequences = file_otolia_values.intersection(file_sen3_values) 

# Find sequences unique to otolia
otolia_only_sequences = file_otolia_values.difference(file_sen3_values)

# Find sequences unique to file2
sen3_only_sequences = file_sen3_values.difference(file_otolia_values)

# Convert the results back to DataFrames if needed
common_df = pd.DataFrame({'sequence': list(common_sequences)})
otolia_only_df = pd.DataFrame({'sequence': list(otolia_only_sequences)})
sen3_only_df = pd.DataFrame({'sequence': list(sen3_only_sequences)})

#print results
print("In total we have:",len(common_df), "kmers common between otolia and Sen3")
print("In total we have:",len(otolia_only_df),"kmers specific to otolia")
print("In total we have:",len(sen3_only_df),"kmers specific to sen3-sephir")

# You can save the DataFrames to files if needed
common_df.to_csv('common_kmers_otolia_sen3.csv', index=False)
otolia_only_df.to_csv('otolia_specific_kmers.csv', index=False)
sen3_only_df.to_csv('Sen3_specific_kmers.csv', index=False)

#Now comparing common_kmers_otolia_sen3 with laura, kuba, to remove the kmers obtained as a result of intersection and only keep remianing (unique) ones in all three datasets. 
file_kuba_values = process_file(file_kuba)
file_laura_values = process_file(file_laura) 

## Use common_df from the previous comparison and convert it into a set
common_df_sen3_otolia=set(common_df['sequence'])  

# Find sequences that are common between laura, kuba, and otolia, and exclude them from the three files plus also the results of intersection for this particular analysis
intersection_sequences_KLC = file_kuba_values.intersection(file_laura_values, common_df_sen3_otolia)

# Remove the intersection sequences from LAURA
file_kuba_values1=pd.DataFrame(file_kuba_values)
file_kuba_values1.columns=["sequence"]
kuba_df = file_kuba_values1[~file_kuba_values1['sequence'].isin(intersection_sequences_KLC)]

# Remove the intersection sequences from KUBA
file_laura_values1=pd.DataFrame(file_laura_values)
file_laura_values1.columns=["sequence"]
laura_df = file_laura_values1[~file_laura_values1['sequence'].isin(intersection_sequences_KLC)]

# Remove the intersection sequences from common_df_SEN3_OTOLIA
common_df_sen3_otolia1=pd.DataFrame(common_df_sen3_otolia)
common_df_sen3_otolia1.columns=["sequence"]
common_df_sen3_otolia_comp2 = common_df_sen3_otolia1[~common_df_sen3_otolia1['sequence'].isin(intersection_sequences_KLC)]

#print results
print("In total we excluded:",len(intersection_sequences_KLC), "common kmers of sen3 and otolia shared with kuba and laura")
print("In total we have:",len(kuba_df),"kmers specific to kuba")
print("In total we have:",len(laura_df),"kmers specific to laura")
print("In total we have:",len(common_df_sen3_otolia_comp2),"kmers common between otlia and sen3-sephir not shared with laura and kuba")
print("In the initial analysis we had",len(common_df), "kmers shared between sen3-sephir and otolia, and in comparison two this number reduced to:",len(common_df_sen3_otolia_comp2),"by removing laura and kuba kmers shared with sen3-otolia list")

#SAVE RESULTS INTO FILES
intersection_sequences_KLC.to_csv("Excluded_kmers.csv",index=False)
kuba_df.to_csv("kuba_only_excluded_interesection_kmers.csv",index=False)
laura_df.to_csv("laura_only_excluded_interesection_kmers.csv",index=False)
common_df_sen3_otolia_comp2.to_csv("Selected_kmers_intersection_sen3_otolia_not_shared_with_kuba_laura.csv",index=False)

