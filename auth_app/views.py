import json

from django.http import HttpResponse
from django.shortcuts import render

from auth_app.models import Users


# Create your views here.
def login_user(request):
    if request.method == 'POST':
        req_body = json.loads(request.body.decode('utf-8'))
        print(req_body)
        try:
            for user in Users.objects:
                if user.email == req_body['email'] and user.password == req_body['password']:
                    return HttpResponse('User Logged In Successfully', status=200)

            return HttpResponse('User Not Found', status=404)

        except Exception as e:
            return HttpResponse("Error \n %s" % e, status=500)


def register_user(request):
    if request.method == 'POST':
        try:
            req_body = json.loads(request.body.decode('utf-8'))
            user = Users()
            user.name = req_body['fullname']
            user.email = req_body['email']
            user.password = req_body['password']

            user.save()
            return HttpResponse('User Registered Successfully')

        except Exception as e:
            return HttpResponse("Error \n %s" % e, status=500)

    else:
        return render(request, '404.html')
