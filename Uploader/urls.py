from django.urls import path
from django.contrib import admin
from . import views

app_name = 'Uploader'
urlpatterns = [
    path('admin/', admin.site.urls),
    path('', views.SDFMultView, name='upload')
]
