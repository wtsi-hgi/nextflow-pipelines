import gzip, io, os, sys
import pandas as pd

fofn_countsfiles = pd.read_csv(sys.argv[1])
lib_counts = sys.argv[2]

for index, row in fofn_countsfiles.iterrows():
    print(row['count_file'])
    if index == 0:
        concatenated = pd.read_csv(row['count_file'], sep='\t')
        concatenated['samplename']=row['samplename']
    else:
        to_concatenate = pd.read_csv(row['count_file'], sep='\t')
        to_concatenate['samplename']=row['samplename']
        concatenated = pd.concat([concatenated, to_concatenate], ignore_index=True)
    print(concatenated.head())
        
count_matrix = concatenated.pivot(index='Guide Sequence', columns='samplename', values='Count')       
count_matrix.to_csv(lib_counts + '.count_matrix.txt', sep='\t')
    

# matrix with gene symbols added
for index, row in fofn_countsfiles.iterrows():
    print(row['count_file'].replace('.counts.txt','.genes.counts.txt'))
    if index == 0:
        concatenated = pd.read_csv(row['count_file'].replace('.counts.txt','.genes.counts.txt'), sep='\t')
        concatenated['samplename']=row['samplename']
    else:
        to_concatenate = pd.read_csv(row['count_file'].replace('.counts.txt','.genes.counts.txt'), sep='\t')
        to_concatenate['samplename']=row['samplename']
        concatenated = pd.concat([concatenated, to_concatenate], ignore_index=True)
    print(concatenated.head())
        
count_matrix = concatenated.pivot_table(index=['Guide Sequence','Gene'], columns='samplename', values='Count').reset_index()       
count_matrix.to_csv(lib_counts + '.genes.count_matrix.txt', sep='\t', index = False)
