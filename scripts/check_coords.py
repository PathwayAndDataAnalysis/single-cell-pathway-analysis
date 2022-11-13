import numpy as np
import pandas as pd


def check_coords(cell_array, umap_coords_path):
    umap_coords = pd.read_table(umap_coords_path, index_col=0)
    traject = umap_coords[umap_coords.index.isin(cell_array)]
    ord_traject = traject.loc[cell_array]
    return(ord_traject)


def get_line_map(trajectory, umap_coords_path):
    umap_coords = pd.read_table(umap_coords_path, index_col=0)
    lines = np.zeros((2, len(trajectory)))
    i = 0
    for t in trajectory:
        lines[0, i] = umap_coords.loc[t][0]
        lines[1, i] = umap_coords.loc[t][1]
        i = i+1

    return lines
