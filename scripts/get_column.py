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


def subset_marker_expression(exp_data_dir, marker_gene_list):
    exp_dict = dict()

    with open(exp_data_dir, "r") as exp_file:

        for e in exp_file:
            key_genes, *lines = e.split('\t')

            for g in marker_gene_list:
                if g in key_genes:
                    exp_dict[key_genes] = lines

    return exp_dict
