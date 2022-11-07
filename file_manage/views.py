import os
import json
from datetime import datetime

from django.http import HttpResponse, JsonResponse


def safe_open_w(path):
    """ Open "path" for writing, creating any parent directories as needed.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return open(path, 'w')


def httpMethodError(expMethod, recMethod):
    return HttpResponse('Expected method: ' + expMethod + " but received method: " + recMethod, status=404)


# Create your views here.
def upload_file(request):
    if request.method == 'POST':
        try:
            data_str = request.body.decode('utf-8')
            # print(data_str)

            # extract filename
            name_data = ((data_str.split("\n", 4)[1]).split(";")[-1]).split('"')[1]
            print("name_data", name_data)

            data_str = data_str.split('\n', 4)[4]
            data_str = data_str[:data_str.rfind('\n')]
            data_str = data_str[:data_str.rfind('\n')]

            # print(data_str)
            filename = "data/kisan@gmail.com/" + name_data
            with safe_open_w(filename) as f:
                f.write(data_str)

            return HttpResponse("File Upload", status=200)

        except Exception as e:
            print(e)
            return HttpResponse('Exception: ' + str(e), status=404)
    else:
        print(request)
        return httpMethodError("POST", request.method)


def convert_date(timestamp):
    d = datetime.utcfromtimestamp(timestamp)
    formatted_date = d.strftime('%d %b %Y')
    return formatted_date


def get_all_files(request):
    if request.method == "GET":
        try:
            files = []
            with os.scandir("data/kisan@gmail.com/") as entries:
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


def delete_file(request):
    if request.method == "POST":
        try:
            req_body = json.loads(request.body.decode('utf-8'))
            print(req_body)

            file_to_delete = "data/kisan@gmail.com/" + req_body['file']
            # If file exists, delete it
            if os.path.isfile(file_to_delete):
                os.remove(file_to_delete)
            else:
                print("Error: %s file not found" % file_to_delete)
            return HttpResponse("Success", status=200)
        except Exception as e:
            print(e)
            return HttpResponse('Exception: ' + str(e), status=404)
    else:
        print(request)
        return httpMethodError("POST", request.method)
