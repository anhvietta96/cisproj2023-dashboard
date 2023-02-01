from django.urls import path
from django.contrib import admin
from . import views

app_name = 'compounds'

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', views.main_compound_view, name='compounds'),
    path('<str:inchi_key>', views.molecule_single_view)
]
