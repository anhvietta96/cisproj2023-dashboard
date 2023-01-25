from django.shortcuts import render
from compounds.models import Molecule,MoleculeSet
import json
from django.http import HttpResponse
from dashboard.settings import MEDIA_ROOT
import os

num_property_list = Molecule.objects.get_num_attr()
display_num_property_list = [property.replace('_',' ').title() for property in num_property_list]
available_chart_types = ['Scatter','Bubble','Histogram']


def ChartOptions(request):
    response = {'data':{}}

    response['data']['property_list']=display_num_property_list

    all_set = MoleculeSet.objects.all()
    set_list=[str(set) for set in all_set]
    response['data']['set_list']=set_list

    response['data']['chart_types']=available_chart_types
    return render(request,'chart_options.html',response)


def ChartResult(request):
    ChartOptions = {'legend':[],'name':[],'data':[],'image':{}}

    property_list = Molecule.objects.get_all_attr()

    post = request.POST
    chart_type = int(post.get('chart-type')[0])

    all_set = MoleculeSet.objects.all()

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
                for mol in q:
                    data_count+=1
                    mol_inf = [getattr(mol, num_property_list[x_axis]),
                               getattr(mol, num_property_list[y_axis])]
                    for property in property_list:
                        if property != num_property_list[x_axis] and property != num_property_list[y_axis] and property not in ['image']:
                            mol_inf.append(getattr(mol,property))
                        if property == 'image':
                            ChartOptions['image'][mol.inchi_key] = mol.image.url
                    data.append(mol_inf)
                ChartOptions['data'].append(data)

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
                    mol_inf = [getattr(mol, num_property_list[x_axis])]
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

    directory = os.path.join(MEDIA_ROOT,'exported_data/csv/')
    if not os.path.exists(directory):
        os.makedirs(directory)
    filepath = os.path.join(directory,filename)

    file = open(filepath,'w+')
    file.writelines(data_to_write)
    file.close()
    with open(filepath,'r') as f:
        file_data = f.read()

    response = HttpResponse(file_data,content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename={}'.format(filename)

    return response
