import dataclasses
from compounds.filters import attr_name_mapping
from django.shortcuts import render
import compounds.filters
from .models import Molecule, MoleculeSet
from django.http import HttpRequest
from django.core.exceptions import ObjectDoesNotExist


@dataclasses.dataclass
class Set:
    set_name: str
    mol_list: list


def main_compound_view(request: HttpRequest):
    set_list = [Set(mol_set.set_name, mol_set.molecules.all())
                for mol_set in MoleculeSet.objects.all()]
    if request.method == "POST":
        lower = request.POST.get("lower")
        upper = request.POST.get("upper")
        pattern = request.POST.get("pattern")
        attr = request.POST.get("attr")
        if attr:
            if lower and upper:
                filter_query = [float(lower), float(upper)]
            elif pattern:
                filter_query=pattern
            else:
                return render(request, 'compounds/compounds.html', {'set_list': set_list})

            model_attr = attr_name_mapping[attr]
            compounds.filters.filter_sets(model_attr, filter_query, set_list)
    return render(request, 'compounds/compounds.html', {'set_list': set_list})


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
