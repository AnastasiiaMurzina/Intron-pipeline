from gtfparse import read_gtf
import pandas as pd
import numpy as np
from copy import deepcopy
import csv

df = read_gtf("gencode.v29.basic.annotation.gtf")

df_reduced = df.loc[lambda x: x['gene_type'] == 'protein_coding']
df_reduced = df.loc[lambda x: x['feature'] == 'exon']

dfr_direct = df_reduced.loc[lambda x: x['strand'] == '+']
dfr_reversed = df_reduced.loc[lambda x: x['strand'] == '-']

dfr_direct = dfr_direct.groupby(['gene_id', 'transcript_id', 'gene_type']).count()
[dfr_direct.pop(i) for i in dfr_direct.columns[1::]]

dfr_direct = dfr_direct.rename(index=str, columns={"seqname": "count_exon"})

idx_d = dfr_direct.groupby(['gene_id'])['count_exon'].transform(max) == dfr_direct['count_exon']
dfr_direct = dfr_direct[idx_d]
dfr_direct = dfr_direct.reset_index()
dfr_direct = dfr_direct.drop_duplicates('gene_id')

exond = df_reduced[df_reduced['transcript_id'].isin(dfr_direct['transcript_id'])]
exond = exond.sort_values(by=['transcript_id'])
dfr_direct = dfr_direct.sort_values(by=['transcript_id'])


flag_prev = ''
noverlaps = []
for i in exond.iterrows():
  if flag_prev=='':
    info = deepcopy(i[1])
    flag_prev = i[1]['transcript_id']
    id_inters = [[i[1]['start'], i[1]['end']]]
  elif flag_prev == i[1]['transcript_id']:
    id_inters.append([i[1]['start'], i[1]['end']])
  else:
    id_inters.sort(key=lambda x: x[0])
    s = [id_inters[0]]
    for inx in id_inters[1::]:
      top = s[-1]
      if top[1]<inx[0]: s.append(inx)
      elif top[1]<inx[1]:
        s[-1][1] = max(top[1], inx[1])
    n_exs = len(s)
    for infs_ix, ins in enumerate(s):
      noverlaps.append(pd.DataFrame.from_dict({'seqname': info['seqname'], 'transcript_id': info['transcript_id'], 'strand': '+', 'feature': 'exon', 'gene_id': info['gene_id'], 'gene_type': info['gene_type'], 'exon_number': infs_ix+1, 'start': ins[0], 'end': ins[1], 'count_exon': n_exs}, orient='index'))
    flag_prev = i[1]['transcript_id']
    info = deepcopy(i[1])
    id_inters = [[i[1]['start'], i[1]['end']]]
id_inters.sort(key=lambda x: x[0])
s = [id_inters[0]]
for inx in id_inters[1::]:
  top = s[-1]
  if top[1]<inx[0]: s.append(inx)
  elif top[1]<inx[1]:
    s[-1][1] = max(top[1], inx[1])
n_exs = len(s)
for infs_ix, ins in enumerate(s):
      noverlaps.append(pd.DataFrame.from_dict({'seqname': info['seqname'], 'transcript_id': info['transcript_id'], 'strand': '+', 'feature': 'exon', 'gene_id': info['gene_id'], 'gene_type': info['gene_type'], 'exon_number': infs_ix+1, 'start': ins[0], 'end': ins[1], 'count_exon': n_exs}, orient='index'))
with open('exonds_noverlapes.tsv', 'w') as csvfile:
  noverlaps[0].T.to_csv(csvfile, index=list(noverlaps[0].index), sep='\t')
  for elem in noverlaps[1::]:
    elem.T.to_csv(csvfile, sep='\t', header=0)


dfr_reversed = dfr_reversed.groupby(['gene_id', 'transcript_id', 'gene_type']).count()
[dfr_reversed.pop(i) for i in dfr_reversed.columns[1::]]

dfr_reversed = dfr_reversed.rename(index=str, columns={"seqname": "count_exon"})

idx_r = dfr_reversed.groupby(['gene_id'])['count_exon'].transform(max) == dfr_reversed['count_exon']
dfr_reversed = dfr_reversed[idx_r]
dfr_reversed = dfr_reversed.reset_index()
dfr_reversed = dfr_reversed.drop_duplicates('gene_id')

exonr = df_reduced[df_reduced['transcript_id'].isin(dfr_reversed['transcript_id'])]
exonr = exonr.sort_values(by=['transcript_id'])
dfr_reversed = dfr_reversed.sort_values(by=['transcript_id'])

flag_prev = ''
rnoverlaps = []
for i in exonr.iterrows():
  if flag_prev=='':
    info = deepcopy(i[1])
    flag_prev = i[1]['transcript_id']
    id_inters = [[i[1]['start'], i[1]['end']]]
  elif flag_prev == i[1]['transcript_id']:
    id_inters.append([i[1]['start'], i[1]['end']])
  else:
    id_inters.sort(key=lambda x: x[0])
    s = [id_inters[0]]
    for inx in id_inters[1::]:
      top = s[-1]
      if top[1]<inx[0]: s.append(inx)
      elif top[1]<inx[1]:
        s[-1][1] = max(top[1], inx[1])
    n_exs = len(s)
    for infs_ix, ins in enumerate(s):
      rnoverlaps.append(pd.DataFrame.from_dict({'seqname': info['seqname'], 'transcript_id': info['transcript_id'], 'strand': '+', 'feature': 'exon', 'gene_id': info['gene_id'], 'gene_type': info['gene_type'], 'exon_number': int(infs_ix)+1, 'start': int(ins[0]), 'end': int(ins[1]), 'count_exon': int(n_exs)}, orient='index'))
    flag_prev = i[1]['transcript_id']
    info = deepcopy(i[1])
    id_inters = [[i[1]['start'], i[1]['end']]]
id_inters.sort(key=lambda x: x[0])
s = [id_inters[0]]
for inx in id_inters[1::]:
  top = s[-1]
  if top[1]<inx[0]: s.append(inx)
  elif top[1]<inx[1]:
    s[-1][1] = max(top[1], inx[1])
n_exs = len(s)
for infs_ix, ins in enumerate(s):
      rnoverlaps.append(pd.DataFrame.from_dict({'seqname': info['seqname'], 'transcript_id': info['transcript_id'], 'strand': '+', 'feature': 'exon', 'gene_id': info['gene_id'], 'gene_type': info['gene_type'], 'exon_number': int(infs_ix)+1, 'start': int(ins[0]), 'end': int(ins[1]), 'count_exon': int(n_exs)}, orient='index'))
with open('exonrs_noverlapes.tsv', 'w') as csvfile:
  rnoverlaps[0].T.to_csv(csvfile, index=list(rnoverlaps[0].index), sep='\t')
  for elem in rnoverlaps[1::]:
    elem.T.to_csv(csvfile, sep='\t', header=0)

###########Nooverlapped exond ans exonr are ready#####################################
ndf = pd.read_csv('exonds_noverlapes.tsv', sep='\t', header=0)
dfd_in = pd.DataFrame(data={'gene_id': ndf['gene_id'].where(ndf['exon_number'] != ndf['count_exon']),
                           'transcript_id': ndf['transcript_id'].where(ndf['exon_number']!= ndf['count_exon']),
                           # 'intron_name': ndf['transcript_id'].where(ndf['exon_number'] != ndf['count_exon'])+'.'+ndf['exon_number'].where(ndf['exon_number'] != ndf['count_exon']),
                           'chr': ndf['seqname'].where(ndf['exon_number']!=ndf['count_exon']), 'strand': ndf['strand'].where(ndf['exon_number']!=ndf['count_exon']),
                           'exon_number': ndf['exon_number'].where(ndf['exon_number']!=ndf['count_exon']),
                           'start': ndf['end'].where(ndf['exon_number']!=ndf['count_exon'])+1, 'end': ndf['start'].where(ndf['exon_number'] != 1).shift(-1)-1,
                           'length': ndf['start'].where(ndf['exon_number'] != 1).shift(-1)-1-ndf['end'].where(ndf['exon_number']!=ndf['count_exon'])})

dfd_in = dfd_in.dropna()
dfd_in.to_csv('intronds_noverlapes.tsv', sep='\t', header=True)

rdf = pd.read_csv('exonrs_noverlapes.tsv', sep='\t', header=0)
rfd_in = pd.DataFrame(data={'gene_id': ndf['gene_id'].where(ndf['exon_number'] != ndf['count_exon']),
                           'transcript_id': ndf['transcript_id'].where(ndf['exon_number']!= ndf['count_exon']),
                           # 'intron_name': ndf['transcript_id'].where(ndf['exon_number'] != ndf['count_exon'])+'.'+ndf['exon_number'].where(ndf['exon_number'] != ndf['count_exon']),
                           'chr': ndf['seqname'].where(ndf['exon_number']!=ndf['count_exon']), 'strand': ndf['strand'].where(ndf['exon_number']!=ndf['count_exon']),
                           'exon_number': ndf['exon_number'].where(ndf['exon_number']!=ndf['count_exon']),
                           'start': ndf['end'].where(ndf['exon_number']!=ndf['count_exon'])+1, 'end': ndf['start'].where(ndf['exon_number'] != 1).shift(-1)-1,
                           'length': ndf['start'].where(ndf['exon_number'] != 1).shift(-1)-1-ndf['end'].where(ndf['exon_number']!=ndf['count_exon'])})

rfd_in = dfd_in.dropna()
rfd_in.to_csv('intronds_noverlapes.tsv', sep='\t', header=True)
#####################Two databases of intron positions##########################################
