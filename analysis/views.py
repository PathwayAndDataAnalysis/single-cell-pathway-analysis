import json
import os
import shutil
import pandas as pd
from django.http import JsonResponse, HttpResponse

from SingleCellPathwayAnalysis.views import httpMethodError, exceptionInRequest
from analysis.models import Analysis, TASK_STATUS
from file_manage.views import convert_date, safe_open_w
from pipeline.views import runPipeline


def get_all_analysis(request):
    """
    Return response listing all the analysis

    @param request:
    @return:
    """
    if request.method == "GET":
        try:
            files = []
            with os.scandir("data/kisan@gmail.com/analysis/") as entries:
                for entry in entries:
                    if entry.is_dir():
                        fileInfo = entry.stat()
                        creation_date = convert_date(fileInfo.st_ctime)

                        sizeUnit = "KB"
                        fileSize = fileInfo.st_size / 1024

                        if fileSize > 1024:
                            fileSize = fileSize / 1024
                            sizeUnit = "MB"

                        if fileSize > 1024:
                            fileSize = fileSize / 1024
                            sizeUnit = "GB"

                        # Open parameter file
                        param_file = 'data/kisan@gmail.com/analysis/' + entry.name + "/" + "params.json"
                        with open(param_file) as f:
                            params = json.load(f)

                        # Get the status of the analysis
                        filteringStatus = Analysis.objects(analysisName=params['analysisName']).first().isFilteringDone
                        pcaStatus = Analysis.objects(analysisName=params['analysisName']).first().isPCADone
                        umapStatus = Analysis.objects(analysisName=params['analysisName']).first().isUMAPDone
                        allStatus = Analysis.objects(analysisName=params['analysisName']).first().isAllDone
                        errorMessage = Analysis.objects(analysisName=params['analysisName']).first().errorMessage

                        files.append({"analysisName": entry.name,
                                      "analysisSize": str(round(fileSize, 3)) + " " + sizeUnit,
                                      "creationDate": creation_date,
                                      "analysisParams": params,

                                      "filteringStatus": filteringStatus,
                                      "pcaStatus": pcaStatus,
                                      "umapStatus": umapStatus,
                                      "allStatus": allStatus,
                                      "errorMessage": errorMessage,
                                      })

            return JsonResponse({'analysis': files})

        except Exception as e:
            return exceptionInRequest(e)

    else:
        return httpMethodError("GET", request.method)


def delete_analysis(request):
    """
    Delete the analysis
    @param request:
    @return:
    """
    if request.method == "POST":
        try:
            req_body = json.loads(request.body.decode('utf-8'))
            print('req_body', req_body)
            file_to_delete = "data/kisan@gmail.com/analysis/" + req_body['analysis']
            shutil.rmtree(file_to_delete, ignore_errors=True, )

            # delete from database
            Analysis.objects(analysisName=req_body['analysis']).delete()

            return HttpResponse("Success", status=200)
        except Exception as e:
            return exceptionInRequest(e)
    else:
        return httpMethodError("POST", request.method)


def create_param_file(req_body):
    """
    Create a new param file params.json file
    @param req_body:
    @return:
    """
    param_file = {
        "analysisName": req_body['analysisName'],
        "dataMatrixFile": req_body['dataMatrixFile'],
        "organism": req_body['organism'],
        "metaDataFile": req_body['metaDataFile'],

        "isFilterCells": req_body['isFilterCells'],
        "isFilterGenes": req_body['isFilterGenes'],
        "isQCFilter": req_body['isQCFilter'],
        "isNormalizeData": req_body['isNormalizeData'],
        "isUseLogTransform": req_body['isUseLogTransform'],

        "usePCA": req_body['usePCA'],
        "useUMAP": req_body['useUMAP'],
        "normalizationScale": req_body['normalizationScale'],
    }

    if req_body['isFilterCells']:
        param_file['minNumOfCells'] = req_body['minNumOfCells']

    if req_body['isFilterGenes']:
        param_file['minNumOfGenes'] = req_body['minNumOfGenes']

    if req_body['isQCFilter']:
        param_file['qcFilterPercent'] = req_body['qcFilterPercent']

    if req_body['isNormalizeData']:
        param_file['normalizationScale'] = req_body['normalizationScale']

    if req_body['usePCA']:
        param_file['pcaCount'] = req_body['pcaCount']
    else:
        param_file['pcaFile'] = req_body['pcaFile']

    if req_body['useUMAP']:
        param_file['n_neighbors'] = req_body['n_neighbors']
        param_file['min_dist'] = req_body['min_dist']
        param_file['metric'] = req_body['metric']
    else:
        param_file['umapFile'] = req_body['umapFile']

    return param_file


def update_database(analysis_name: str,
                    filtering_status: TASK_STATUS = TASK_STATUS.DEFAULT.value,
                    pca_status: TASK_STATUS = TASK_STATUS.DEFAULT.value,
                    umap_status: TASK_STATUS = TASK_STATUS.DEFAULT.value,
                    completed_status: TASK_STATUS = TASK_STATUS.DEFAULT.value, ):
    """
    Update the database with the status of the analysis
    @param analysis_name:
    @param filtering_status:
    @param pca_status:
    @param umap_status:
    @param completed_status:
    """
    # Create a new entry in the database
    analysis = Analysis()
    analysis.analysisName = analysis_name
    analysis.isFilteringDone = filtering_status
    analysis.isPCADone = pca_status
    analysis.isUMAPDone = umap_status
    analysis.isAllDone = completed_status
    analysis.errorMessage = ""
    analysis.save()


def run_analysis(request):
    """
    Run the analysis
    @param request:
    @return:
    """
    if request.method == "POST":
        try:
            req_body = json.loads(request.body.decode('utf-8'))
            print('req_body', req_body)

            # Create params.json file content
            new_param_file = create_param_file(req_body)

            # Convert into json
            param_file_json = json.dumps(new_param_file, indent=4)

            analysis_file = 'data/kisan@gmail.com/analysis/' + req_body['analysisName'] + '/params.json'
            with safe_open_w(analysis_file) as f:
                f.write(param_file_json)

            # update database
            update_database(req_body['analysisName'])

            # run pipeline
            file_path = 'data/kisan@gmail.com/files' + '/'
            analysis_path = 'data/kisan@gmail.com/analysis/' + req_body['analysisName'] + '/'
            new_param_file['exp_path'] = file_path + new_param_file['dataMatrixFile']
            new_param_file['metadata_path'] = file_path + new_param_file['metaDataFile']
            new_param_file['out_path'] = analysis_path + "_filtered.tsv"

            new_param_file['pca_exp_path'] = analysis_path + "_filtered.tsv"
            new_param_file['pca_out_path'] = analysis_path + "_pca.tsv"

            new_param_file['umap_exp_path'] = analysis_path + "_pca.tsv"
            new_param_file['umap_out_path'] = analysis_path + "_umap.tsv"

            new_param_file['umap_clustered_path'] = analysis_path + "_umap_clustered.tsv"
            runPipeline(new_param_file)

            return HttpResponse(status=200)

        except Exception as e:
            return exceptionInRequest(e)

    else:
        return httpMethodError("GET", request.method)


def update_analysis(request):
    """
    Update the analysis
    @param request:
    @return:
    """
    if request.method == "POST":
        try:
            req_body = json.loads(request.body.decode('utf-8'))
            print('req_body: ', req_body)

            analysis_name = req_body['analysisName']
            old_analysis_name = req_body['oldAnalysisName']

            # Create params.json file content
            new_param_file = create_param_file(req_body)

            # Convert into json
            param_file_json = json.dumps(new_param_file, indent=4)

            # Check if analysis name is changed, if is not changed then update the same file
            if analysis_name == old_analysis_name:
                path = 'data/kisan@gmail.com/analysis/' + old_analysis_name + '/params.json'
                with safe_open_w(path) as f:
                    f.write(param_file_json)
            else:
                # Rename analysis folder
                old_path = 'data/kisan@gmail.com/analysis/' + old_analysis_name
                new_path = 'data/kisan@gmail.com/analysis/' + analysis_name
                os.rename(old_path, new_path)

                new_path = new_path + '/params.json'

                new_param_file['analysisName'] = analysis_name
                param_file_json = json.dumps(new_param_file, indent=4)
                with safe_open_w(new_path) as f:
                    f.write(param_file_json)

            # update database
            update_database(new_param_file['analysisName'])

            # run pipeline
            file_path = 'data/kisan@gmail.com/files' + '/'
            analysis_path = 'data/kisan@gmail.com/analysis/' + req_body['analysisName'] + '/'
            new_param_file['exp_path'] = file_path + new_param_file['dataMatrixFile']
            new_param_file['out_path'] = analysis_path + "_filtered.tsv"

            new_param_file['pca_exp_path'] = analysis_path + "_filtered.tsv"
            new_param_file['pca_out_path'] = analysis_path + "_pca.tsv"

            new_param_file['umap_exp_path'] = analysis_path + "_pca.tsv"
            new_param_file['umap_out_path'] = analysis_path + "_umap.tsv"
            runPipeline(new_param_file)

            return HttpResponse(status=200)

        except Exception as e:
            return exceptionInRequest(e)

    else:
        return httpMethodError("GET", request.method)


def get_coordinates(request):
    """
    Get the coordinates of the analysis
    @param request:
    @return:
    """
    if request.method == "POST":
        try:
            req_body = json.loads(request.body.decode('utf-8'))
            print('req_body', req_body)

            analysis_name = req_body['analysisName']

            keys = ["Cell", "UMAP1", "UMAP2", "ClusterID"]
            data_to_return = []
            coordinate_file = "data/kisan@gmail.com/analysis/" + analysis_name + "/_umap_clustered.tsv"
            with open(coordinate_file, 'r') as f:
                lines = f.readlines()[1:]
                for line in lines:
                    ll = [i.strip() for i in line.split('\t')]
                    dd = {}
                    for i, l in enumerate(ll):
                        dd[keys[i]] = l
                    data_to_return.append(dd)
            return JsonResponse({'data': data_to_return}, status=200)

        except Exception as e:
            return exceptionInRequest(e)

    else:
        return httpMethodError("POST", request.method)


def get_metadata_columns(request):
    """
    Get the metadata columns
    @param request:
    @return:
    """
    if request.method == "POST":
        try:
            req_body = json.loads(request.body.decode('utf-8'))
            analysis_name = req_body['analysisName']

            meta_file = "data/kisan@gmail.com/analysis/" + analysis_name + "/params.json"
            with open(meta_file, 'r') as f:
                meta_data = json.load(f)
                meta_file = meta_data['metaDataFile']

                mt = pd.read_csv("data/kisan@gmail.com/files/" + meta_file, sep='\t')
                cols = list(mt.columns)[1:]
                return JsonResponse({'data': cols}, status=200)

        except Exception as e:
            return exceptionInRequest(e)
    else:
        return httpMethodError("POST", request.method)


def get_data_with_metadata_columns(request):
    """
    Get the data with metadata columns
    @param request:
    @return:
    """
    if request.method == "POST":
        try:
            req_body = json.loads(request.body.decode('utf-8'))

            analysis_name = req_body['analysisName']
            column_name = req_body['columnName']

            param_file = "data/kisan@gmail.com/analysis/" + analysis_name + "/params.json"
            umap_path = "data/kisan@gmail.com/analysis/" + analysis_name + "/_umap.tsv"
            out_path = "data/kisan@gmail.com/analysis/" + analysis_name + "/_umap_clustered.tsv"

            with open(param_file, 'r') as f:
                meta_data = json.load(f)
                meta_file = meta_data['metaDataFile']

                # Read and manage metadata file
                mt = pd.read_csv("data/kisan@gmail.com/files/" + meta_file, sep='\t')
                hd = list(mt.columns)
                hd[0] = 'cell_id'
                mt.columns = hd
                # replace - with .
                mt['cell_id'] = mt['cell_id'].str.replace('-', '.')

                # Read and manage umap
                coord = pd.read_csv(umap_path, sep='\t')
                hd = list(coord.columns)
                hd[0] = 'cell_id'
                coord.columns = hd

                # merge metadata and umap
                result = coord.merge(mt, how='left', on='cell_id').drop_duplicates()
                result = result[['cell_id', 'umap comp. 1', 'umap comp. 2', column_name]]
                result.to_csv(out_path, sep='\t', index=False)

                # Response data
                # keys = ["Cell", "UMAP1", "UMAP2", "ClusterID"]
                # data_to_return = []
                # coordinate_file = "data/kisan@gmail.com/analysis/" + analysis_name + "/_umap_clustered.tsv"
                # with open(coordinate_file, 'r') as f:
                #     lines = f.readlines()[1:]
                #     for line in lines:
                #         ll = [i.strip() for i in line.split('\t')]
                #         dd = {}
                #         for i, l in enumerate(ll):
                #             dd[keys[i]] = l
                #         data_to_return.append(dd)
                #
                return JsonResponse({'data': "success"}, status=200)

        except Exception as e:
            return exceptionInRequest(e)
    else:
        return httpMethodError("POST", request.method)
