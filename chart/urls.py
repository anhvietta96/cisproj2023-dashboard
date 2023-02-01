from django.urls import path
from .views import ChartOptions, ChartResult, Export_CSV

urlpatterns = [
    path('', ChartOptions, name='chart_options'),
    path('result/', ChartResult, name='chart_result'),
    path('result/tsv', Export_TSV, name='tsv_download'),
]
