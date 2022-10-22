from django.http import HttpResponse
from django.shortcuts import render

from auth_app.models import Users


# Create your views here.
def login_user(request):
    return render(request, 'login.html', {'name': 'Kisan'})


def register_user(request):
    user = Users.objects.create(
        name="request.POST['name']",
        email="request.POST['email']",
        password="request.POST['password']"
        # name=request.POST['name'],
        # email=request.POST['email'],
        # password=request.POST['password']
    )
    user.save()
    return HttpResponse("User Registered!!!")
