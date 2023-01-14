from django.shortcuts import render
from django.views.generic.edit import CreateView
from compounds.models import Molecule,MoleculeSet
'''from .models import PropertyChoice
from .forms import PropertyForm'''
import json
from django.http import HttpResponse
from dashboard.settings import MEDIA_ROOT
import os

num_property_list = Molecule.objects.get_num_attr()
display_num_property_list = [property.replace('_',' ').title() for property in num_property_list]
available_chart_types = ['Scatter','Bubble','Histogram']

class Chart:
    pass

# Create your views here.
def ChartOptions(request):
    response = {'data':{}}
    
    response['data']['property_list']=display_num_property_list

    all_set = MoleculeSet.objects.all()
    set_list=[str(set) for set in all_set]
    response['data']['set_list']=set_list
    
    response['data']['chart_types']=available_chart_types
    return render(request,'chart_options.html',response)
'''
def ChartResult(request):
    
    x_axis = int(request.POST.get('x-axis')[0])-1
    y_axis = int(request.POST.get('y-axis')[0])-1
    q = Molecule.objects.all()
    
    properties = Molecule.__dict__["__doc__"]
    property_list = properties[9:].replace(',','').replace(')','').split()

    ChartOptions = {"chart":{"height":350,"type":"scatter","zoom":{"enabled":True,"type":"xy"}},
    "xaxis": {
    "tickAmount": 10,"labels" : {"show":True}, "title":{"text":property_list[x_axis]}
    },"yaxis": {
    "tickAmount": 7,"labels" : {"show":True}, "title":{"text":property_list[y_axis]}}}

    data = []
    info = {"name":property_list[x_axis]+"/"+property_list[y_axis]}
    info["data"] = []
    for mol in q:
        info["data"].append([getattr(mol,property_list[x_axis]),getattr(mol,property_list[y_axis])])
    data.append(info)
    
    ChartOptions["series"]=data
    result = {'title':'Created Chart','options':ChartOptions}
    return render(request,'chart_result.html',result)
'''

def ChartResult(request):
    ChartOptions = {'legend':[],'name':[],'data':[],'image':[]}

    property_list = Molecule.objects.get_all_attr()

    post = request.POST
    print(post)
    chart_type = int(post.get('chart-type')[0])

    all_set = MoleculeSet.objects.all()
    set_list=[str(set) for set in all_set]

    x_axis = int(post.get('x-axis')[0])-1
    
    post_dict = post.dict()
    set_name_list = []

    lexicographic_position = []

    data_count = 0

    if(chart_type == 1 or chart_type == 2):
        y_axis = int(post.get('y-axis')[0])-1

        for key in post_dict.keys():
            if key.startswith('set_'):
                data = []
                set_num = int(post.get(key)[0])-1
                set = all_set[set_num]
                set_name_list.append(set.set_name)
                q = set.molecules.all()
                img_list = []
                for mol in q:
                    data_count+=1
                    mol_inf = []
                    mol_inf.append(getattr(mol,num_property_list[x_axis]))
                    mol_inf.append(getattr(mol,num_property_list[y_axis]))
                    for property in property_list:
                        if property != num_property_list[x_axis] and property != num_property_list[y_axis] and property not in ['image']:
                            mol_inf.append(getattr(mol,property))
                        if property == 'image':
                            img_list.append(mol.image.url)
                    data.append(mol_inf)
                ChartOptions['data'].append(data)
                ChartOptions['image'].append(img_list)

        ChartOptions['header']=[num_property_list[x_axis],num_property_list[y_axis]]
        for property in property_list:
            if property != num_property_list[x_axis] and property != num_property_list[y_axis] and property not in ['image']:
                ChartOptions['header'].append(property)
        for i,property in enumerate(ChartOptions['header']):
            if property in Molecule.objects.get_str_attr():
                lexicographic_position.append(i)
        ChartOptions['header']=[header.replace('_',' ').title() for header in ChartOptions['header']]
        ChartOptions['legend']=set_name_list
        ChartOptions['lexicographic_position'] = lexicographic_position
        ChartOptions['size'] = data_count
        title = display_num_property_list[x_axis]+"/"+display_num_property_list[y_axis]+" Chart for Group "+", ".join(set_name_list)
    elif chart_type == 3:
        for key in post_dict.keys():
            if key.startswith('set_'):
                data = []
                set_num = int(post.get(key)[0])-1
                set = all_set[set_num]
                set_name_list.append(set.set_name)
                q = set.molecules.all()
                for mol in q:
                    data_count+=1
                    mol_inf = []
                    mol_inf.append(getattr(mol,num_property_list[x_axis]))
                    for property in property_list:
                        if property != num_property_list[x_axis] and property not in ['image']:
                            mol_inf.append(getattr(mol,property))
                    data.append(mol_inf)
                ChartOptions['data'].append(data)

        ChartOptions['header']=[num_property_list[x_axis]]
        for property in property_list:
            if property != num_property_list[x_axis] and property not in ['image']:
                ChartOptions['header'].append(property)
        for i,property in enumerate(ChartOptions['header']):
            if property in Molecule.objects.get_str_attr():
                lexicographic_position.append(i)

        ChartOptions['header']=[header.replace('_',' ').title() for header in ChartOptions['header']]
        ChartOptions['legend']=set_name_list
        ChartOptions['lexicographic_position'] = lexicographic_position
        ChartOptions['size'] = data_count
        title = display_num_property_list[x_axis]+" Histogram for Group "+", ".join(set_name_list)
    
    ChartOptions['type']=chart_type
    result = {'title':title,'options':ChartOptions}
    return render(request,'chart_result.html',result)

def Export_CSV(request):
    inchikey_collection = json.loads(request.POST['export-csv-val'])

    filename = 'exported_csv_{}.csv'.format(request.POST['csrfmiddlewaretoken'])
    
    property_list = Molecule.objects.get_all_export_attr()
    data_to_write = ['\t'.join(property_list)+'\n']
    for inchikey in inchikey_collection:
        data_string = ""
        mol = Molecule.objects.filter(inchi_key__exact=inchikey)
        for property in property_list:
            data_string += str(getattr(mol[0],property))
            data_string += '\t'
        data_string = data_string[:-1] + '\n'
        data_to_write.append(data_string)
    print(data_to_write)
    filepath = os.path.join(MEDIA_ROOT,'exported_data/csv/',filename)
    file = open(filepath,'w+')
    file.writelines(data_to_write)
    file.close()
    with open(filepath,'r') as f:
        file_data = f.read()

    response = HttpResponse(file_data,content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename={}'.format(filename)

    return response

def Export_SDF(request):
    return