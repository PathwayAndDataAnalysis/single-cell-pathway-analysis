import threading
import time
from django.http import JsonResponse, HttpResponse

from pipeline.models import ThreadTask


# Create your views here.

def runPipeline(request):
    runPCA(request)
    runUMAP(request)


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
    time.sleep(10)
    pca_task = ThreadTask.objects.get(pk=id)
    pca_task.is_done = True
    pca_task.save()
    print("Finished PCA Process id", id)


# UMAP Pipeline
def runUMAP(request):
    print("Starting UMAP Process")
    umap_task = ThreadTask.objects.create(task_name="umap_process")
    umap_task.save()
    umap_thread = threading.Thread(target=umapRunningTask, args=[umap_task.id], daemon=True)
    umap_thread.start()
    return JsonResponse({'id': umap_task.id})


def checkUMAP(request, id):
    umap_task = ThreadTask.objects.get(pk=id)
    return JsonResponse({'is_done': umap_task.is_done})


def umapRunningTask(id):
    print("Starting UMAP Process id", id)
    time.sleep(10)
    umap_task = ThreadTask.objects.get(pk=id)
    umap_task.is_done = True
    umap_task.save()
    print("Finished UMAP Process id", id)





# def startThreadTask(request):
#     task = ThreadTask.objects.create(task_name="test")
#     task.save()
#     thread1 = threading.Thread(target=longRunningTask, args=[task.id], daemon=True)
#     thread1.start()
#     return JsonResponse({'id': task.id})
#
#
# def checkThreadTask(request, id):
#     task = ThreadTask.objects.get(pk=id)
#     return JsonResponse({'is_done': task.is_done})
#
#
# def longRunningTask(id):
#     print("Starting long running task", id)
#     time.sleep(10)
#     task = ThreadTask.objects.get(pk=id)
#     task.is_done = True
#     task.save()
#     print("Finished long running task", id)



