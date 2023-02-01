from django.shortcuts import render

import compounds.filters
from .models import Molecule, MoleculeSet
from django.http import HttpRequest
from django.core.exceptions import ObjectDoesNotExist


def main_compound_view(request: HttpRequest):
    sets = MoleculeSet.objects.all()
    return render(request, 'compounds/compounds.html', {'sets': sets})

attr_name_mapping = {
    "Formula": "molecular_formula",
    "logP": "log_p",
    "Weight": "molecular_weight",
    "H-acceptors": "num_h_acceptors",
    "H-donors": "num_h_donors",
    "Rotatable bonds": "num_rotatable_bonds",
    "Rings": "num_rings"
}

def filter_compound_view(request: HttpRequest):
    lower = request.POST.get("lower")
    upper = request.POST.get("upper")
    pattern = request.POST.get("pattern")
    attr = request.POST.get("attr")
    if lower != None and upper != None:
        filter_query = [int(lower), int(upper)]
    elif pattern != None:
        filter_query=pattern 

    sets = MoleculeSet.objects.all()
    if attr != None:
        model_attr = attr_name_mapping[attr]
        sets = compounds.filters.filter_sets(model_attr, filter_query)
    return render(request, 'compounds/filtered.html', {'sets': sets})

def molecule_single_view(request: HttpRequest, inchi_key: str):
    try:
        molecule = Molecule.objects.get(pk=inchi_key)
    except ObjectDoesNotExist:
        return render(request, 'compounds/compounds.html')

    return render(request, 'compounds/compounds.html',
                  {'molecules': [molecule]})

def search_results(request):
    query = request.GET.get("search_query")
    q = Molecule.objects.filter(name__icontains=query)
    search_results = []
    property_list = Molecule.objects.get_all_attr()

    property_list.remove('image')
    property_list = ['image'] + property_list
    display_property_list = [property.replace('_',' ').title() for property in property_list]

    for mol in q:
        value_list = []
        for property in property_list:
            if property != 'image':
                value_list.append(getattr(mol, property))
            else:
                value_list.append(mol.image.url)
        search_results.append(value_list)

    data = {'header': display_property_list, 'table': search_results}
    return render(request, 'search.html', data)
