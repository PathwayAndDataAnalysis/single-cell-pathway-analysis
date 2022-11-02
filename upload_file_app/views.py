import os

from django.http import HttpResponse
from django.shortcuts import render
import csv

from upload_file_app.models import UserData


# Create your views here.
def upload_file(request):
    if request.method == 'POST':
        try:
            user_data = UserData()
            user_data.user_email = "kisan@gmail.com"
            user_data.cell_data = str(request.body)

            user_data.save()

            return HttpResponse("File Upload", status=200)

        except Exception as e:
            print(e)
            return HttpResponse('Exception: ' + str(e), status=404)
    else:
        print(request)
        return HttpResponse('Not Found', status=404)
