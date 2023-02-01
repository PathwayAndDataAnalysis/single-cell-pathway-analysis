import pandas as pd

def subset_marker_expression(exp_data_dir, marker_gene_list):

  df = pd.read_csv(exp_data_dir, sep='\t', index_col=0)

  df = df.loc[marker_gene_list]
  df = df.sum(axis=0)
  return df


data_dir = './_filtered.tsv'
# gene_list = ['Tcea1', 'Sema4c', 'Il18r1',
#              'Mfsd6', 'Hibch', '9430016H08Rik', 'Fam126b', 'Creb1', 'Rqcd1', 'Zfp142', 'Bcs1l', 'Rnf25', 'Ttll4', 'Wnt10a', 'Cnppd1',
#              'Rhobtb2', 'Rhbdd1', 'Mff', 'Agfg1', 'Trip12', 'Fbxo36', 'Sp140', 'Cab39', 'Itm2c', 'Psmd1', 'Ncl', 'Ptma', 'Pde6d']

gene_list = ['Tcea1', 'Sema4c', 'Mfsd6', 'Hibch', '9430016H08Rik']

print("Testing subset_marker_expression()...")

val = subset_marker_expression(data_dir, gene_list)

df = pd.DataFrame(val)
df.to_csv('test_output.tsv', sep='\t', index=True, header=False)

print(val)
