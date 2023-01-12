import threading

import pandas as pd
import scanpy as sc
import umap
from django.http import JsonResponse
from pandas import DataFrame

from analysis.models import Analysis
from pipeline.models import ThreadTask


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
    adata = sc.read(norm_exp_path).T

    # Find highly variable genes and subtract them
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata = adata[:, adata.var.highly_variable]

    # Scale then perform PCA analysis
    sc.pp.scale(adata)
    sc.tl.pca(adata, svd_solver='arpack')

    # Write the k PCA components in to a file
    pc_table = adata.obsm['X_pca'][:, 0:n_pcs]
    pd.DataFrame(data=pc_table, index=adata.obs_names).to_csv(out_path, sep="\t")


# input path should be the output of sc_PCA function
def run_umap(norm_exp_path: str, input_path: str, output_path: str, metric: str = 'euclidean', min_dist: float = 0.1,
             n_neighbors: int = 15) -> None:
    """
    This function gets PCA components as an input file and find distances using those components.
    It outputs file containing 2D locations using UMAP function.
    @param norm_exp_path:
    @param input_path:
    @param output_path:
    @param metric:
    @param min_dist:
    @param n_neighbors:
    """
    adata = sc.read(norm_exp_path).T
    # Read the file containing k PCA components
    input_data = pd.read_table(input_path, index_col=0)
    # transpose is needed because pca output is transposed
    input_t = input_data

    # Run UMAP
    reducer = umap.UMAP(metric=metric, min_dist=min_dist, n_neighbors=n_neighbors, random_state=0)
    embedding = reducer.fit_transform(input_t)

    # Write the output of UMAP as a file
    umap_df: DataFrame = pd.DataFrame(data=embedding, index=adata.obs_names, columns=['umap comp. 1', 'umap comp. 2'])
    umap_df.to_csv(output_path, sep="\t")


# Add Clustering
def addClustering(metadata_path, umap_path, out_path):
    # Read and manage metadata
    metadata = pd.read_csv(metadata_path, sep='\t')
    hd = list(metadata.columns)
    hd[0] = 'cell_id'
    metadata.columns = hd
    # replace - with .
    metadata['cell_id'] = metadata['cell_id'].str.replace('-', '.')

    # Read and manage umap
    coord = pd.read_csv(umap_path, sep='\t')
    hd = list(coord.columns)
    hd[0] = 'cell_id'
    coord.columns = hd

    # merge metadata and umap
    result = coord.merge(metadata, how='left', on='cell_id').drop_duplicates()
    result = result[['cell_id', 'umap comp. 1', 'umap comp. 2', 'seurat_clusters']]
    result.to_csv(out_path, sep='\t', index=False)


# Create your views here.
def runPipeline(analysis_info):
    runTask(analysis_info)


# Pipeline
def runTask(analysis_info):
    print("Starting PCA Process")
    task = ThreadTask.objects.create(task_name="pca_process")
    task.save()
    thread = threading.Thread(target=processRunningTask, args=[task.id, analysis_info], daemon=True)
    thread.start()
    return JsonResponse({'id': task.id})


def checkTaskThread(request, id):
    task = ThreadTask.objects.get(pk=id)
    return JsonResponse({'is_done': task.is_done})


def processRunningTask(id, analysis_info):
    print("Starting PCA Process id", id)

    # Do Filtering
    print('analysis_info', analysis_info)
    filter_and_normalize(
        species=analysis_info['organism'],
        apply_min_cells=analysis_info['isFilterCells'],
        apply_min_genes=analysis_info['isFilterGenes'],
        apply_mt_filter=analysis_info['isQCFilter'],
        apply_normalization=analysis_info['isNormalizeData'],
        apply_log=analysis_info['isUseLogTransform'],
        exp_path=analysis_info['exp_path'],
        out_path=analysis_info['out_path'],
        min_genes=int(analysis_info['minNumOfGenes']),
        min_cells=int(analysis_info['minNumOfCells']),
        target_sum=int(analysis_info['normalizationScale']),
        percent_mt=float(analysis_info['qcFilterPercent']),
    )
    Analysis.objects(analysisName=analysis_info['analysisName']).update(isFilteringDone=True)
    print("Filtering Done")

    # Do PCA
    run_pca(
        norm_exp_path=analysis_info['pca_exp_path'],
        out_path=analysis_info['pca_out_path'],
        n_pcs=int(analysis_info['pcaCount']),
    )
    Analysis.objects(analysisName=analysis_info['analysisName']).update(isPCADone=True)
    print("PCA Done")

    # Do UMAP
    run_umap(
        norm_exp_path=analysis_info['pca_exp_path'],
        input_path=analysis_info['umap_exp_path'],
        output_path=analysis_info['umap_out_path'],
        metric=analysis_info['metric'],
        min_dist=float(analysis_info['min_dist']),
        n_neighbors=int(analysis_info['n_neighbors']),
    )
    Analysis.objects(analysisName=analysis_info['analysisName']).update(isUMAPDone=True)
    print("UMAP Done")

    # Do Clustering
    addClustering(
        metadata_path=analysis_info['metadata_path'],
        umap_path=analysis_info['umap_out_path'],
        out_path=analysis_info['umap_clustered_path'],
    )

    # Update Database
    Analysis.objects(analysisName=analysis_info['analysisName']).update(isAllDone=True)

    task = ThreadTask.objects.get(pk=id)
    task.is_done = True
    task.save()
    print("Everything done, Process id", id)
