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

def run_umap(norm_exp_path: str, input_path: str, out_path: str) -> None:
    """

    @param norm_exp_path:
    @param input_path:
    @param out_path:
    @param n_pcs:
    """
    adata = sc.read(norm_exp_path)

    # Read the file containing k PCA components
    input_data = pd.read_table(input_path, index_col=0)
    adata.obsm['X_pca'] = input_data

    # Compute nearest neighbors
    sc.pp.neighbors(adata)
    sc.tl.umap(adata, n_components=np.shape(adata.obsm['X_pca'])[1])

    #
    umap_table = adata.obsm['X_umap']
    pd.DataFrame(data=umap_table).to_csv(out_path)