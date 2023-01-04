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
    time.sleep(5)
    print("Filtering Done")
    # Update Database


    # Do PCA
    time.sleep(5)
    print("PCA Done")

    # Update Database

    # Do UMAP
    time.sleep(5)
    print("UMAP Done")
    # Update Database
    analysis = Analysis.objects(analysis_name="test").update(upsert=True, set__isDone=True)
    # analysis.is_everything_done = True
    # analysis.save()

    pca_task = ThreadTask.objects.get(pk=id)
    pca_task.is_done = True
    pca_task.save()
    print("Finished PCA Process id", id)

# UMAP Pipeline
# def runUMAP(request):
#     print("Starting UMAP Process")
#     umap_task = ThreadTask.objects.create(task_name="umap_process")
#     umap_task.save()
#     umap_thread = threading.Thread(target=umapRunningTask, args=[umap_task.id], daemon=True)
#     umap_thread.start()
#     return JsonResponse({'id': umap_task.id})
#
#
# def checkUMAP(request, id):
#     umap_task = ThreadTask.objects.get(pk=id)
#     return JsonResponse({'is_done': umap_task.is_done})
#
#
# def umapRunningTask(id):
#     print("Starting UMAP Process id", id)
#     time.sleep(10)
#     umap_task = ThreadTask.objects.get(pk=id)
#     umap_task.is_done = True
#     umap_task.save()
#     print("Finished UMAP Process id", id)
