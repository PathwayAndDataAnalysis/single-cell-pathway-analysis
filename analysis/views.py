import os

from django.http import JsonResponse, HttpResponse

from file_manage.views import httpMethodError, convert_date


def get_all_analysis(request):
    if request.method == "GET":
        try:
            files = []
            with os.scandir("data/kisan@gmail.com/analysis/") as entries:
                for entry in entries:
                    print(entry)
                    fileInfo = entry.stat()
                    uploadDate = convert_date(fileInfo.st_ctime)
                    fileSize = fileInfo.st_size / 1024

                    print("fileInfo", fileInfo)
                    files.append({"fileName": entry.name, "fileSize": round(fileSize, 3), "uploadDate": uploadDate})

            return JsonResponse({'files': files})

        except Exception as e:
            print(e)
            return HttpResponse('Exception: ' + str(e), status=404)

    else:
        print(request)
        return httpMethodError("GET", request.method)
