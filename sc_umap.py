import pandas as pd
import umap

# input path should be the output of sc_PCA function
def sc_umap (input_path, output_path):
    input_data = pd.read_table(input_path,index_col=0)
    input_T = input_data.T
    um = umap.UMAP()
    X_fit = um.fit(input_T)
    X_umap = um.transform(input_T)
    umap_df = pd.DataFrame(data = X_umap, columns = ['umap comp. 1', 'umap comp. 2'])
    umap_df.to_csv(output_path,sep="\t")
    return (umap_df)