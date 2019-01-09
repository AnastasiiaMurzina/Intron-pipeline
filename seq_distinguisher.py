from Bio import SeqIO
import pandas as pd

comlement = {'A': 'T', 'C':'G', 'T':'A','G':'C', 'N':'N'}
sequences = SeqIO.parse('GRCh38.primary_assembly.genome.fa', 'fasta')
dfd_in = pd.read_csv('intronds_noverlapes.tsv', sep='\t', header=0)
dfr_in = pd.read_csv('intronrs_noverlapes.tsv', sep='\t', header=0)

# with open('GRCh38_introns.fasta', 'w') as f:                                                                                                                                                               
with open('G_introns.fasta', 'w') as f:                                                                                                                                                               
    for i in sequences:
        for j in dfd_in[dfd_in['chr'] == i.id].iterrows():
            f.write('>' +str(j[1]['transcript_id'])+'.'+str(int(j[1]['length']))+'.'+str(int(j[1]['exon_number']))+'\n')
            ll = int(j[1]['length'])
            for k in range(ll//60):
                f.write(str(i.seq[int(j[1]['start'])-1+60*k:min(len(i), int(j[1]['start'])-1+60*(k+1)):]))
                f.write('\n')
            f.write('\n')
        for j in dfr_in[dfr_in['chr'] == i.id].iterrows():
            f.write('>' +str(j[1]['transcript_id']+'.'+str(j[1]['length'])) +'\n')
            for k, ix in enumerate(str(i.seq[int(j[1]['start'])-1:int(j[1]['end'])][::-1]):
                    f.write(comlement[ix])
                    if k != 0 and k % 59 == 0:
                        f.write('\n')
