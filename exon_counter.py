from gtfparse import read_gtf
import pandas as pd
from copy import deepcopy
import sys

def filters(df, **kwags):
  df_reduced = df
  for i, j in kwags.items():
    df_reduced = df_reduced.loc[lambda x: x[i] == j]
  return df_reduced

def merge_overlaped(intervals):
  intervals.sort(key=lambda x: x[0])
  s = [intervals[0]]
  for inx in intervals[1::]:
    top = s[-1]
    if top[1]<inx[0]: s.append(inx)
    elif top[1]<inx[1]:
      s[-1][1] = max(top[1], inx[1])
  return s

def nooverlaped(df, exon_file, strand='+'):
  df_cur = deepcopy(df)
  df_cur = df_cur.groupby(['gene_id', 'transcript_id', 'gene_type']).count()
  [df_cur.pop(i) for i in df_cur.columns[1::]]
  df_cur = df_cur.rename(index=str, columns={"seqname": "count_exon"})
  df_cur = df_cur[df_cur.groupby(['gene_id'])['count_exon'].transform(max) == df_cur['count_exon']].reset_index().drop_duplicates('gene_id').sort_values(by=['transcript_id'])
  exons = df[df['transcript_id'].isin(df_cur['transcript_id'])].sort_values(by=['transcript_id'])
  flag_prev = ''
  noverlaps = []
  info = deepcopy(exons.iloc[0])
  flag_prev = exons.iloc[0]['transcript_id']
  id_inters = [[exons.iloc[0]['start'], exons.iloc[0]['end']]]
  for i in exons[1::].iterrows():
    if flag_prev == i[1]['transcript_id']:
      id_inters.append([i[1]['start'], i[1]['end']])
    else:
      s = merge_overlaped(id_inters)      
      n_exs = len(s)
      for infs_ix, ins in enumerate(s):
        noverlaps.append(pd.DataFrame.from_dict({'seqname': info['seqname'], 'transcript_id': info['transcript_id'], 'strand': strand, 'feature': 'exon', 'gene_id': info['gene_id'], 'gene_type': info['gene_type'], 'exon_number': infs_ix+1, 'start': ins[0], 'end': ins[1], 'count_exon': n_exs}, orient='index'))
      flag_prev = i[1]['transcript_id']
      info = deepcopy(i[1])
      id_inters = [[i[1]['start'], i[1]['end']]]
  s = merge_overlaped(id_inters)
  n_exs = len(s)
  for infs_ix, ins in enumerate(s):
        noverlaps.append(pd.DataFrame.from_dict({'seqname': info['seqname'], 'transcript_id': info['transcript_id'], 'strand': strand, 'feature': 'exon', 'gene_id': info['gene_id'], 'gene_type': info['gene_type'], 'exon_number': infs_ix+1, 'start': ins[0], 'end': ins[1], 'count_exon': n_exs}, orient='index'))
  with open(exon_file, 'w') as csvfile:
    noverlaps[0].T.to_csv(csvfile, index=list(noverlaps[0].index), sep='\t')
    for elem in noverlaps[1::]:
      elem.T.to_csv(csvfile, sep='\t', header=0)

def exons_to_introns(fname, ofname):
    ndf = pd.read_csv(fname, sep='\t', header=0)
    dfd_in = pd.DataFrame(data={'gene_id': ndf['gene_id'].where(ndf['exon_number'] != ndf['count_exon']),
                             'transcript_id': ndf['transcript_id'].where(ndf['exon_number']!= ndf['count_exon']),
                             'chr': ndf['seqname'].where(ndf['exon_number']!=ndf['count_exon']), 'strand': ndf['strand'].where(ndf['exon_number']!=ndf['count_exon']),
                             'exon_number': ndf['exon_number'].where(ndf['exon_number']!=ndf['count_exon']),
                             'start': ndf['end'].where(ndf['exon_number']!=ndf['count_exon'])+1, 'end': ndf['start'].where(ndf['exon_number'] != 1).shift(-1)-1,
                             'length': ndf['start'].where(ndf['exon_number'] != 1).shift(-1)-1-ndf['end'].where(ndf['exon_number']!=ndf['count_exon'])})

    dfd_in = dfd_in.dropna()
    dfd_in.to_csv(ofname, sep='\t', header=True)
  


if __name__ == '__main__':
  gtf_annotation = sys.argv[1]  
  df = read_gtf(gtf_annotation)
  i = 2
  fils = {'feature': 'exon', 'gene_type': 'protein_coding'}
  while i < len(sys.argv) and '=' in sys.argv[i]:
    what, which = sys.argv[i].split('=')
    fils.update({what: which})
    i+=1

  dfile = gtf_annotation[:min(6, len(gtf_annotation))] if len(sys.argv)==i else sys.argv[i]+'d_noverlaped.tsv'
  rfile = gtf_annotation[:min(6, len(gtf_annotation))] if len(sys.argv)==i else sys.argv[i]+'r_noverlaped.tsv'
  
  df_reduced = filters(df, **fils)

  dfr_direct = df_reduced.loc[lambda x: x['strand'] == '+']
  dfr_reversed = df_reduced.loc[lambda x: x['strand'] == '-']

  nooverlaped(dfr_direct, 'exonds_noverlapes.tsv')
  nooverlaped(dfr_reversed, 'exonrs_noverlapes.tsv', strand='-')

  exons_to_introns('exonds_noverlapes.tsv', dfile)
  exons_to_introns('exonrs_noverlapes.tsv', rfile)
