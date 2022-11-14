from django.urls import path

from . import views

urlpatterns = [
    path('get_all/', views.get_all_analysis, name='get_all_analysis'),
    path('run/', views.run_analysis, name='run_analysis'),
    path('delete/', views.delete_analysis, name='delete_analysis'),
]
