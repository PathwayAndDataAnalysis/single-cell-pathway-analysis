import pandas as pd
import scanpy as sc


def filter_and_normalize(species: str, apply_min_cells: bool, apply_min_genes: bool, apply_mt_filter: bool,
                         apply_normalization: bool, apply_log: bool, exp_path: str, out_path: str,
                         min_genes: int = 200, min_cells: int = 10, target_sum: int = 1e4,
                         percent_mt: float = 10) -> None:
    adata = sc.read(exp_path).T

    if apply_min_cells:
        sc.pp.filter_cells(adata, min_genes=min_genes)  # 200
    if apply_min_genes:
        sc.pp.filter_genes(adata, min_cells=min_cells)  # 3

    if apply_mt_filter:
        if species == "Mouse":
            adata.var['mt'] = adata.var_names.str.startswith('mt-')
        elif species == "Human":
            adata.var['mt'] = adata.var_names.str.startswith('MT-')
        else:
            raise Exception("This species is not supported: " + species)

        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        adata = adata[adata.obs.pct_counts_mt < percent_mt, :]

    if apply_normalization:
        sc.pp.normalize_total(adata, target_sum=target_sum)
    if apply_log:
        sc.pp.log1p(adata)

    pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names).to_csv(out_path)
