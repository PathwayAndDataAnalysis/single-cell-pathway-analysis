import json
import os
import shutil

from django.http import JsonResponse, HttpResponse

from analysis.models import Analysis
from file_manage.views import httpMethodError, convert_date, safe_open_w
from pipeline.views import runPipeline


def get_all_analysis(request):
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

                        files.append({"analysisName": entry.name,
                                      "analysisSize": str(round(fileSize, 3)) + " " + sizeUnit,
                                      "creationDate": creation_date,
                                      "analysisParams": params,

                                      "filteringStatus": filteringStatus,
                                      "pcaStatus": pcaStatus,
                                      "umapStatus": umapStatus,
                                      "allStatus": allStatus,
                                      })

            return JsonResponse({'analysis': files})

        except Exception as e:
            print(e)
            return HttpResponse('Exception: ' + str(e), status=404)

    else:
        print(request)
        return httpMethodError("GET", request.method)


def delete_analysis(request):
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
            print(e)
            return HttpResponse('Exception: ' + str(e), status=404)
    else:
        print(request)
        return httpMethodError("POST", request.method)


def create_param_file(req_body):
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


def run_analysis(request):
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
            analysis = Analysis()
            analysis.analysisName = new_param_file['analysisName']
            analysis.isFilteringDone = False
            analysis.isPCADone = False
            analysis.isUMAPDone = False
            analysis.isAllDone = False
            analysis.save()

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
            print(e)
            return HttpResponse('Exception: ' + str(e), status=404)

    else:
        print(request)
        return httpMethodError("GET", request.method)


def update_analysis(request):
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
            # analysis = Analysis()
            # analysis.analysisName = new_param_file['analysisName']
            # analysis.isFilteringDone = False
            # analysis.isPCADone = False
            # analysis.isUMAPDone = False
            # analysis.isAllDone = False
            # analysis.save()

            # run pipeline
            # runPipeline(request)

            return HttpResponse(status=200)

        except Exception as e:
            print(e)
            return HttpResponse('Exception: ' + str(e), status=404)

    else:
        print(request)
        return httpMethodError("GET", request.method)


def get_coordinates(request):
    if request.method == "POST":
        try:
            req_body = json.loads(request.body.decode('utf-8'))
            print('req_body', req_body)

            analysis_name = req_body['analysisName']

            keys = ["Cell", "UMAP1", "UMAP2"]
            data_to_return = []
            with open("data/kisan@gmail.com/analysis/" + analysis_name + "/_umap.tsv", 'r') as f:
                lines = f.readlines()[1:]
                for line in lines:
                    ll = [i.strip() for i in line.split('\t')]
                    dd = {}
                    for i, l in enumerate(ll):
                        dd[keys[i]] = l
                    data_to_return.append(dd)

            return JsonResponse({'data': data_to_return}, status=200)

        except Exception as e:
            print(e)
            return HttpResponse('Exception: ' + str(e), status=404)

    else:
        print(request)
        return httpMethodError("POST", request.method)
