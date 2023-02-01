from typing import Dict, Any
import pandas as pd


def get_meta_column(meta_data_dir: str, column: str) -> None:
    """

    @param meta_data_dir:
    @param column:
    @return:
    """
    meta_data = pd.read_table(meta_data_dir, index_col=0)
    meta_column_dict = dict(zip(meta_data[column].index, meta_data[column].values))
    return meta_column_dict


# def subset_marker_expression(exp_data_dir, marker_gene_list):
#     exp_dict = dict()
#
#     with open(exp_data_dir, "r") as exp_file:
#
#         for e in exp_file:
#             key_genes, *lines = e.split('\t')
#
#             for g in marker_gene_list:
#                 if g in key_genes:
#                     exp_dict[key_genes] = lines
#
#     return exp_dict

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
