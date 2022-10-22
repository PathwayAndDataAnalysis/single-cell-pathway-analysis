from django.shortcuts import render
from django.http import HttpResponse


# Create your views here.
def login_user(request):
    return render(request, 'login.html', {'name': 'Kisan'})
