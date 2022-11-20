from django.shortcuts import render

def main_compound_view(request):
    return render(request, 'compounds/compounds.html')
