from django.shortcuts import render
from django.template.loader import get_template


def index(request):
    return render(request, 'index.html')


def searchQuery(request):
    data = {}
    if request.method == 'GET':
        search_str = request.GET.get('search_str')
        data = {'data': search_str}
        template = get_template('landing.html')
    return render(request, 'landing.html', context=data)
