from Bio import SeqIO
import pandas as pd
import sys

comlement = {'A': 'T', 'C':'G', 'T':'A','G':'C', 'N':'N'}
sequences = SeqIO.parse(sys.argv[1], 'fasta')

dfd_in = pd.read_csv(sys.argv[2], sep='\t', header=0)
dfr_in = pd.read_csv(sys.argv[3], sep='\t', header=0)

with open(sys.argv[4], 'w') as f:                                                                                                                                                               
    for i in sequences:
        for j in dfd_in[dfd_in['seqname'] == i.id].iterrows():
            ll = int(j[1]['end']-j[1]['start'])
            f.write('>' +str(j[1]['transcript_id'])+'.'+str(ll)+'.'+str(int(j[1]['exon_number']))+'\n')
            for k in range(ll//60):
                f.write(str(i.seq[int(j[1]['start'])-1+60*k:min(len(i), int(j[1]['start'])-1+60*(k+1)):]))
                f.write('\n')
            f.write('\n')
        for j in dfr_in[dfr_in['seqname'] == i.id].iterrows():
            ll = int(j[1]['end']-j[1]['start'])
            f.write('>' +str(j[1]['transcript_id'])+'.'+str(ll) +'\n')
            for k, ix in enumerate(str(i.seq[int(j[1]['start'])-1:int(j[1]['end']):-1])):
                    f.write(comlement[ix])
                    if k != 0 and k % 59 == 0:
                        f.write('\n')
