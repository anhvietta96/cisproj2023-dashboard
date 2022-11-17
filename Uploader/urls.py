from django.urls import include,path
from django.contrib import admin
from django.conf.urls.static import static
from django.conf import settings
from . import views

urlpatterns = [path('admin/',admin.site.urls),path('',views.upload)]