from django.urls import path
from . import views

urlpatterns = [
    path('upload/', views.upload_file, name='file_upload'),
    path('get_all/', views.get_all_files, name='get_all_files'),
    path('delete/', views.delete_file, name="delete_file")
]
