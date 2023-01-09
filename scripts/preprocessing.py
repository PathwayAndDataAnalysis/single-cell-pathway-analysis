import pandas as pd
import scanpy as sc
import umap
from pandas import DataFrame


def filter_and_normalize(species: str, apply_min_cells: bool, apply_min_genes: bool, apply_mt_filter: bool,
                         apply_normalization: bool, apply_log: bool, exp_path: str, out_path: str,
                         min_genes: int = 200, min_cells: int = 10, target_sum: int = 1e4,
                         percent_mt: float = 10) -> None:
    """
    This function applies several preprocessing steps on the expression matrix. According to user preferences it can
    filter genes by minimum number of expression in cells, filter cells by minimum number of genes, filter cells by
    maximum mitochondrial content, it can scal and log normalize.

    @param species:
    @param apply_min_cells:
    @param apply_min_genes:
    @param apply_mt_filter:
    @param apply_normalization:
    @param apply_log:
    @param exp_path: Expression file path
    @param out_path:
    @param min_genes:
    @param min_cells:
    @param target_sum:
    @param percent_mt:
    """
    adata = sc.read(exp_path).T

    # Filter genes and cells
    if apply_min_cells:
        sc.pp.filter_cells(adata, min_genes=min_genes)  # 200
    if apply_min_genes:
        sc.pp.filter_genes(adata, min_cells=min_cells)  # 3

    # Find mitochondrial gene percentage based on the species, filter them using the given cutoff
    if apply_mt_filter:
        if species == "Mouse":
            adata.var['mt'] = adata.var_names.str.startswith('mt-')
        elif species == "Human":
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
        else:
            raise Exception("This species is not supported: " + species)

        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        adata = adata[adata.obs.pct_counts_mt < percent_mt, :]

    # Apply scaling and log normalization
    if apply_normalization:
        sc.pp.normalize_total(adata, target_sum=target_sum)
    if apply_log:
        sc.pp.log1p(adata)

    # Write the processed expression matrix
    pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names).T.to_csv(out_path, sep="\t")


def run_pca(norm_exp_path: str, out_path: str, n_pcs: int = 10) -> None:
    """
    This function gets normalized expression matrix as an input. First it finds highly variable genes and substracts
    them. Then it performs scaling, and PCA.

    @param norm_exp_path:
    @param out_path:
    @param n_pcs:
    """
    adata = sc.read(norm_exp_path)

    # Find highly variable genes and subtract them
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]

    # Scale then perform PCA analysis
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')

    # Write the k PCA components in to a file
    pc_table = adata.obsm['X_pca'][:, 0:n_pcs]
    pd.DataFrame(data=pc_table).T.to_csv(out_path, sep="\t")


# input path should be the output of sc_PCA function
def run_umap(input_path: str, output_path: str, metric: str = 'euclidean', min_dist: float = 0.1,
             n_neighbors: int = 15) -> None:
    """
    This function gets PCA components as an input file and find distances using those components.
    It outputs file containing 2D locations using UMAP function.

    @param input_path:
    @param output_path:
    @param metric:
    @param min_dist:
    @param n_neighbors:
    """
    # Read the file containing k PCA components
    input_data = pd.read_table(input_path, index_col=0)
    # transpose is needed because pca output is transposed
    input_t = input_data.T

    # Run UMAP
    reducer = umap.UMAP(metric=metric, min_dist=min_dist, n_neighbors=n_neighbors)
    embedding = reducer.fit_transform(input_t)

    # Write the output of UMAP as a file
    umap_df: DataFrame = pd.DataFrame(data=embedding, columns=['umap comp. 1', 'umap comp. 2'])
    umap_df.to_csv(output_path, sep="\t")
