import threading

import pandas as pd
import scanpy as sc
from umap.umap_ import UMAP
from django.http import JsonResponse

from analysis.models import Analysis, TASK_STATUS
from pipeline.models import ThreadTask


def filter_and_normalize(species: str, apply_min_cells: bool, apply_min_genes: bool, apply_mt_filter: bool,
                         apply_normalization: bool, apply_log: bool, exp_path: str, out_path: str,
                         min_genes: int = 200, min_cells: int = 10, target_sum: int = 1e4,
                         percent_mt: float = 10):
    """
    This function applies several preprocessing steps on the expression matrix. According to user preferences it can
    filter genes by minimum number of expression in cells, filter cells by minimum number of genes, filter cells by
    maximum mitochondrial content, it can scal and log normalize.
    @param species: Species of the expression matrix
    @param apply_min_cells: If True, filter cells by minimum number of genes
    @param apply_min_genes: If True, filter genes by minimum number of expression in cells
    @param apply_mt_filter:
    @param apply_normalization:
    @param apply_log:
    @param exp_path: Expression file path
    @param out_path: Output file path for the filtered expression matrix file
    @param min_genes:
    @param min_cells:
    @param target_sum:
    @param percent_mt:
    """
    data = sc.read(exp_path).T

    # Filter genes and cells
    if apply_min_cells:
        sc.pp.filter_cells(data, min_genes=min_genes)  # 200
    if apply_min_genes:
        sc.pp.filter_genes(data, min_cells=min_cells)  # 3

    # Find mitochondrial gene percentage based on the species, filter them using the given cutoff
    if apply_mt_filter:
        if species == "Mouse":
            data.var['mt'] = data.var_names.str.startswith('mt-')
        elif species == "Human":
            data.var['mt'] = data.var_names.str.startswith('MT-')
        else:
            raise Exception("This species is not supported: " + species)

        sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        data = data[data.obs.pct_counts_mt < percent_mt, :]

    # Apply scaling and log normalization
    if apply_normalization:
        sc.pp.normalize_total(data, target_sum=target_sum)
    if apply_log:
        sc.pp.log1p(data)

    # Write the processed expression matrix
    pd.DataFrame(data=data.X, index=data.obs_names, columns=data.var_names.str.upper()).to_csv(out_path, sep="\t")


def run_pca(norm_exp_path: str, out_path: str, n_pcs: int = 10) -> None:
    """
    This function gets normalized expression matrix as an input. First it finds highly variable genes and substracts
    them. Then it performs scaling, and PCA.
    @param norm_exp_path: Normalized/Filtered expression matrix file path
    @param out_path: Output file path for the PCA components
    @param n_pcs: Number of PCA components
    """
    data = sc.read(norm_exp_path)

    # Find highly variable genes and subtract them
    sc.pp.highly_variable_genes(data, min_mean=0.0125, max_mean=3, min_disp=0.5)
    data = data[:, data.var.highly_variable]

    # Scale then perform PCA analysis
    sc.pp.scale(data)
    sc.tl.pca(data, svd_solver='arpack')

    # Write the k PCA components in to a file
    pc_table = data.obsm['X_pca'][:, 0:n_pcs]
    pd.DataFrame(data=pc_table, index=data.obs_names).to_csv(out_path, sep="\t", index_label="CellID")


# input path should be the output of sc_PCA function
def run_umap(pca_path: str, output_path: str, metric: str = 'euclidean', min_dist: float = 0.1,
             n_neighbors: int = 15) -> None:
    """
    This function gets PCA components as an input file and find distances using those components.
    It outputs file containing 2D locations using UMAP function.
    @param pca_path: PCA components file path
    @param output_path: Output file path for the UMAP components
    @param metric:
    @param min_dist:
    @param n_neighbors:
    """
    # Read the file containing k PCA components
    pca_data = pd.read_table(pca_path, index_col=0)

    # Run UMAP
    reducer = UMAP(metric=metric, min_dist=min_dist, n_neighbors=n_neighbors, random_state=0)
    embedding = reducer.fit_transform(pca_data)

    # Write the output of UMAP as a file
    pd.DataFrame(data=embedding, index=pca_data.index, columns=['UMAP1', 'UMAP2']).to_csv(output_path, sep="\t",
                                                                                          index_label="CellID")


# Add Clustering
def addClustering(metadata_path, umap_path, out_path):
    """
    This function gets metadata and UMAP components as an input file and merge them.
    @param metadata_path: Metadata file path
    @param umap_path: UMAP components file path
    @param out_path: Output file path for the merged file
    @return:
    """
    metadata = pd.read_csv(metadata_path, sep='\t', index_col=0)
    metadata.index = metadata.index.str.replace('-', '.')
    metadata.index = metadata.index.str.replace(' ', '.')

    coord = pd.read_csv(umap_path, sep='\t', index_col=0)
    # merge metadata and umap
    result = coord.merge(metadata["seurat_clusters"], how="left", left_index=True, right_index=True)
    result.to_csv(out_path, sep='\t')


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


def killTaskThread(id):
    task = ThreadTask.objects.get(pk=id)
    task.is_done = True
    task.save()
    return JsonResponse({'is_done': task.is_done})


def processRunningTask(id, analysis_info):
    print("Starting PCA Process id", id)
    print('analysis_info', analysis_info)

    current_ana = Analysis.objects(analysisName=analysis_info['analysisName'])
    current_ana.update(isAllDone=TASK_STATUS.RUNNING.value)

    # Do Filtering
    try:
        Analysis.objects(analysisName=analysis_info['analysisName']).update(isFilteringDone=TASK_STATUS.RUNNING.value)
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
        current_ana.update(isFilteringDone=TASK_STATUS.COMPLETED.value)
        print("Filtering Done")
    except Exception as e:
        print("Error in Filtering", e)
        current_ana.update(isFilteringDone=TASK_STATUS.FAILED.value)
        current_ana.update(isAllDone=TASK_STATUS.FAILED.value)
        current_ana.update(errorMessage="Error in Filtering: " + str(e))
        killTaskThread(id)
        return

    # Do PCA
    try:
        Analysis.objects(analysisName=analysis_info['analysisName']).update(isPCADone=TASK_STATUS.RUNNING.value)
        run_pca(
            norm_exp_path=analysis_info['pca_exp_path'],
            out_path=analysis_info['pca_out_path'],
            n_pcs=int(analysis_info['pcaCount']),
        )
        current_ana.update(isPCADone=TASK_STATUS.COMPLETED.value)
        print("PCA Done")
    except Exception as e:
        print("Error in PCA", e)
        current_ana.update(isPCADone=TASK_STATUS.FAILED.value)
        current_ana.update(isAllDone=TASK_STATUS.FAILED.value)
        current_ana.update(errorMessage="Error in PCA: " + str(e))
        killTaskThread(id)
        return

    # Do UMAP
    try:
        Analysis.objects(analysisName=analysis_info['analysisName']).update(isUMAPDone=TASK_STATUS.RUNNING.value)
        run_umap(
            pca_path=analysis_info['umap_exp_path'],
            output_path=analysis_info['umap_out_path'],
            metric=analysis_info['metric'],
            min_dist=float(analysis_info['min_dist']),
            n_neighbors=int(analysis_info['n_neighbors']),
        )
        current_ana.update(isUMAPDone=TASK_STATUS.COMPLETED.value)
        print("UMAP Done")
    except Exception as e:
        print("Error in UMAP", e)
        current_ana.update(isUMAPDone=TASK_STATUS.FAILED.value)
        current_ana.update(isAllDone=TASK_STATUS.FAILED.value)
        current_ana.update(errorMessage="Error in UMAP: " + str(e))
        killTaskThread(id)
        return

    # Do Clustering
    try:
        addClustering(
            metadata_path=analysis_info['metadata_path'],
            umap_path=analysis_info['umap_out_path'],
            out_path=analysis_info['umap_clustered_path'],
        )
        print("Clustering Done")
    except Exception as e:
        print("Error in Clustering", e)
        current_ana.update(errorMessage="Error in Clustering: " + str(e))

    # Update Database
    current_ana.update(isAllDone=TASK_STATUS.COMPLETED.value)

    task = ThreadTask.objects.get(pk=id)
    task.is_done = True
    task.save()
    print("Everything done, Process id", id)
