from django.urls import path,include
from django.contrib import admin
from . import views
from rest_framework import routers

app_name = 'compounds'

router = routers.DefaultRouter()
router.register(r"compounds",views.CompoundViewSet)

urlpatterns = [
    path('admin/', admin.site.urls),
    path('', views.main_compound_view, name='compounds'),
    path('data/',include(router.urls),name='data'),
    path('search/',views.SearchResultsView.as_view(),name='search_url'),
    path('<str:inchi_key>', views.molecule_single_view)
]
