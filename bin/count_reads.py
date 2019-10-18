import gzip, io, os, sys
from Bio import SeqIO
import pandas as pd


if len(sys.argv) not in [4,5,6]:
    print('Usage: count_reads.py gzfile library(pilot/main3/main5/test/June) outfile <(opt)plasmid_only> <(opt)includeG>')
else:
    gzfile = sys.argv[1]
    library = sys.argv[2]     #e.g. rasa_in_plasmid.txt for pilot
    outfile = sys.argv[3]
    plasmid_only = eval(sys.argv[4]) if len(sys.argv) == 5 else False 
    includeG = eval(sys.argv[5]) if len(sys.argv) == 6 else False 

    if library == 'pilot':
        guidefile = 'rasa_in_plasmid.txt' if plasmid_only else 'rasa_guides_pilot.txt'
    elif library == 'main3':
        guidefile = 'main_in_plasmid_3.txt' if plasmid_only else 'main_library_guides_3p.txt'
    elif library == 'main5':
        guidefile = 'main_in_plasmid_5.txt' if plasmid_only else 'main_library_guides_5p.txt'
    elif library == 'June3':
        guidefile = './felicity-other/iPSC_double_guide/June/' + 'June_in_plasmid_Post_3.txt' # mod gn5 Pre to Post
    elif library == 'June5':
        guidefile = './felicity-other/iPSC_double_guide/June/' + 'June_in_plasmid_Post_5.txt' # mod gn5 Pre to Post
    else:
        raise Exception('Unknown library:' + library)

    guides = pd.read_csv(guidefile,sep='\t')
    getSeq = lambda x: x[1:] if len(x) == 20 else x
    guide_counts = {getSeq(x):0 for x in guides['Guide Sequence']}
    guide_type = {getSeq(x):y for (x,y) in zip(guides['Guide Sequence'], guides['Type'])}

    count = 0

    if '.gz' in gzfile: f = gzip.open(gzfile, 'rt') # mod gn5
    else: f = open(gzfile,'rt')

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
    f.close()

    fout = io.open(outfile,'w')
    fout.write(u'Guide Sequence\tCount\tType\n')
    fout.write(u'unmapped\t%d\t\n' % (total-mapped))
    for guide in guide_counts:
        fout.write('%s\t%s\t%s\n' % (guide, guide_counts[guide], guide_type[guide]))
    fout.close()

    fout2 = io.open(outfile + '.pcmapped','w') # mod gn5 
    fout2.write(str(outfile) + ' ' + str(mapped) + ' mapped of ' + str(total) + ' ' + str(mapped*100.0/total if total > 0 else 0.0))
    fout2.close()

    print(outfile, mapped, ' mapped of ', total, ' ', mapped*100.0/total if total > 0 else 0.0)
