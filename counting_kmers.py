#counting total kmers
import pandas as pd
import os
import numpy as np

#reading multiple files
file_path=["OTOLIA_inherited_from_Laura.txt", "Otolia_inherited_from_None.txt"]
dfs=[]
for i in file_paths:
    df=pd.read_csv(i,delimiter="\t")
    dfs.append(df)
#concatenate dfs into 1
combined_df=pd.concat(dfs,ignore_index=True)
chr_counts=combined_df.groupby("chromosome")["unique_kmers"].sum()
print(chr_counts)
