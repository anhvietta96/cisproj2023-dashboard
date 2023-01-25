"""dashboard URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import path, include
from django.conf import settings
from django.conf.urls.static import static
from .views import index
from compounds.views import search_results

urlpatterns = [
      path("", index, name="home"),
      path("admin/", admin.site.urls),
      path("upload/", include('Uploader.urls'), name='upload'),
      path("compounds/", include('compounds.urls'), name='compounds'),
      path("chart/", include('chart.urls'), name='chart'),
      path("search/", search_results, name='search_url'),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
