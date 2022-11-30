from django.urls import path,include
from django.contrib import admin
from .views import ChartOptions,ChartResult

urlpatterns = [
    path('', ChartOptions, name='chart'),
    path('result/',ChartResult,name='chart_result')
]
