from django.urls import path,include
from django.contrib import admin
from .views import ChartOptions,ChartResult,Export_CSV,Export_SDF

urlpatterns = [
    path('', ChartOptions, name='chart_options'),
    path('result/',ChartResult,name='chart_result'),
    path('result/csv',Export_CSV,name='csv_download'),
    path('result/sdf',Export_SDF,name='sdf_download'),
]
