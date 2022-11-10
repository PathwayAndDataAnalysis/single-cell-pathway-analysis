import pandas as pd

def check_coords(cell_array,umap_coords_path):
    umap_coords = pd.read_table(umap_coords_path,index_col=0)
    traject = umap_coords[umap_coords.index.isin(cell_array)]
    return(traject)