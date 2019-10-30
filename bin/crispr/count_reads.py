import gzip, io, os, sys
from Bio import SeqIO
import pandas as pd

flog = io.open('record_log.txt','w')
fastq_file = sys.argv[1]
guide = pd.read_csv(sys.argv[2])
outfile = sys.argv[3]
outfile2 = sys.argv[4]
includeG = eval(sys.argv[5])

flog.write("fastq file: " + str(fastq_file) + str("\n"))
flog.write("guide: " + str(sys.argv[2]) + str("\n"))
flog.write("outfile: " + str(outfile) + str("\n"))
flog.write("outfile2: " + str(outfile2) + str("\n"))
flog.write("includeG: " + str(includeG) + str("\n"))

getSeq = lambda x: x[0:19] if len(x) == 20 else x
#getSeq = lambda x: x[1:] if len(x) == 20 else x
guide_counts = {getSeq(x):0 for x in guide['Guide Sequence']}
guide_gene = {getSeq(x):y for (x,y) in zip(guide['Guide Sequence'], guide['Gene'])}
# guide_type = {getSeq(x):y for (x,y) in zip(guide['Guide Sequence'], guides['Type'])}

if '.gz' in fastq_file: f = gzip.open(fastq_file, 'rt') 
else: f = open(fastq_file,'rt')

count = 0
total, mapped = 0, 0
for record in SeqIO.parse(f, 'fastq'):
    if includeG: seq = str(record.seq)[1:20]
    else: seq = str(record.seq)[:19]
    if seq in guide_counts:
        guide_counts[seq] += 1
        mapped += 1
    else: 
        count += 1
        #if count > 1000 and count < 1100:
        #    print(record.seq)
    total += 1
    if total % 100000 == 0:
        # print("total parsed: " + str(total) + " \n")
        flog.write(str(total) + "\n")
flog.close()
f.close()

fout = io.open(outfile,'w')
fout_gene = io.open(outfile.replace('.counts.txt','.genes.counts.txt'),'w')

fout.write(u'Guide Sequence\tCount\n')
fout_gene.write(u'Guide Sequence\tCount\tGene\n')
#fout_gene.write(u'unmapped\t%d\t\n' % (total-mapped))

for guide in guide_counts:
    fout.write('%s\t%s\n' % (guide, guide_counts[guide]))
    fout_gene.write('%s\t%s\t%s\n' % (guide, guide_counts[guide], guide_gene[guide]))
fout.close()
fout_gene.close()

#pd.concat([
#    pd.read_csv(outfile, sep='\t')],
#    guide[['Guide Sequence','Gene']],
#          on='Guide Sequence', how='left').to_csv(
#              outfile.replace('.counts.txt','.genes.counts.txt'), sep='\t', index=False)


fout2 = io.open(outfile2,'w') # mod gn5 
fout2.write(str(outfile) + ' ' + str(mapped) + ' mapped of ' + str(total) + ' ' + str(mapped*100.0/total if total > 0 else 0.0))
fout2.close()
    
# print(outfile, mapped, ' mapped of ', total, ' ', mapped*100.0/total if total > 0 else 0.0)
print(mapped*100.0/total)
