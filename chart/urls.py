from django.urls import path
from .views import ChartOptions, ChartResult, Export_CSV

urlpatterns = [
    path('', ChartOptions, name='chart_options'),
    path('result/', ChartResult, name='chart_result'),
    path('result/csv', Export_CSV, name='csv_download'),
]
