from django.http import HttpResponse
from django.shortcuts import render

from auth_app.models import Users


# Create your views here.
def login_user(request):
    if request.method == 'POST':
        try:
            for user in Users.objects:
                print(user)
                if user.email == request.headers['Email'] and user.password == request.headers['Password']:
                    return HttpResponse('User Logged In Successfully')

            return HttpResponse('User Not Found')

        except Exception as e:
            return HttpResponse("Error \n %s" % e, status=500)


def register_user(request):
    if request.method == 'POST':
        try:
            user = Users()
            user.name = request.body['name']
            user.email = request.body['Email']
            user.password = request.body['Password']

            user.save()
            return HttpResponse('User Registered Successfully')

        except Exception as e:
            return HttpResponse("Error \n %s" % e, status=500)

    else:
        return render(request, '404.html')
