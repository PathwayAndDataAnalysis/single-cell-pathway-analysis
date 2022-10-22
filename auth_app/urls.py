from django.urls import path

from . import views

# URLConf
urlpatterns = [
    path('', views.login_user, name='login'),
]
