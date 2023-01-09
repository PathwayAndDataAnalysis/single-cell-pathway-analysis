import threading
import time

from django.http import JsonResponse

from analysis.models import Analysis
from pipeline.models import ThreadTask


# Create your views here.

def runPipeline(request):
    runPCA(request)
    # runUMAP(request)


# PCA Pipeline
def runPCA(request):
    print("Starting PCA Process")
    pca_task = ThreadTask.objects.create(task_name="pca_process")
    pca_task.save()
    pca_thread = threading.Thread(target=pcaRunningTask, args=[pca_task.id], daemon=True)
    pca_thread.start()
    return JsonResponse({'id': pca_task.id})


def checkPCAThread(request, id):
    pca_task = ThreadTask.objects.get(pk=id)
    return JsonResponse({'is_done': pca_task.is_done})


def pcaRunningTask(id):
    print("Starting PCA Process id", id)

    # Do Filtering
    time.sleep(10)
    # a = Analysis.objects(analysisName="test")
    # Analysis.objects(analysisName="test").update(isFilteringDone=True)
    print("Filtering Done")

    # Do PCA
    time.sleep(10)
    # Analysis.objects(analysisName="test").update(isPCADone=True)
    print("PCA Done")

    # Do UMAP
    time.sleep(10)
    # Analysis.objects(analysisName="test").update(isUMAPDone=True)
    print("UMAP Done")

    # Update Database
    Analysis.objects(analysisName="test").update(isDone=True)

    pca_task = ThreadTask.objects.get(pk=id)
    pca_task.is_done = True
    pca_task.save()
    print("Finished PCA Process id", id)
