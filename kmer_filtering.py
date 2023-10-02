import pandas as pd

#extracting kmers with high frequency
sen3=pd.read_csv("WU18-3664-HR.WU18-3664-NR.subtract.Saphir.intersect.Deodara.subtract.sgenome.subtract.dump",sep="\t")
sen3.columns=["sequence","Freq"]

final_list=pd.read_csv("Mapped_kmers_info.csv",sep="\t")
final_list.columns=["Chrom","Position","sequence"]

#getting the frequency value for final_list and then sorting them in largest to smallest for each chromosme
final_li=pd.merge(sen3,final_list, on=["sequence"])
dataframes=final_li.groupby("Chromosome",group_keys=False).apply(lambda x:x.sort_values("Position",ascending=True))
print(dataframes["Frequency"].mean())

#print("Sorted merged wu and final calls",len(final_li))
dataframes.to_csv("kmers_cmp.csv",sep="\t",index=False)
