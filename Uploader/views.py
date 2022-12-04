from django.shortcuts import render
from .forms import SDFileForm,SDFileMult
from rdkit import Chem
from dashboard.settings import MEDIA_ROOT
from compounds.models import Molecule

def SDFile_Upload(request):
    form=SDFileForm()
    data={'form':form}
    if request.method == 'POST':
        form=SDFileForm(request.POST,request.FILES)
        if form.is_valid():
            form.save()
            
            if file.endswith('.sdf'):
                files=request.FILES.getlist('document')
                file_path = MEDIA_ROOT + '/uploaded_data/'
                with Chem.SDMolSupplier(file_path+file) as suppl:
                    for substance in suppl:
                        if substance is None:
                            continue
                            subst = Molecule()
                            subst.inchi_key = substance.GetProp('PUBCHEM_IUPAC_INCHIKEY')
                            subst.log_p = substance.GetProp('PUBCHEM_XLOGP3')
                            subst.num_h_acceptors = substance.GetProp('PUBCHEM_CACTVS_HBOND_ACCEPTOR')
                            subst.num_h_donors = substance.GetProp('PUBCHEM_CACTVS_HBOND_DONOR')
                            subst.molecular_mass = substance.GetProp('PUBCHEM_MOLECULAR_WEIGHT')
                            subst.save()
                data['form']=form
                data['existing_uploads']={'headers':['File','Set','Description']}
                data['existing_uploads']['data'] = [[file,1,form['description'].value()]]
            else:
                err_msg='Not a SDFile'
                data['err_msg']=err_msg
    return render(request,'upload.html',data)

def save_file(file):
    with open(MEDIA_ROOT + '/uploaded_data/' + file.name,'wb+') as destination:
        for chunk in file.chunks():
            destination.write(chunk)

def SDFMultView(request):
    form=SDFileMult()
    data={'form':form}
    if request.method == 'POST':
        form=SDFileMult(request.POST,request.FILES)
        if form.is_valid():
            data['existing_uploads']={'headers':['File','Set','Description'],'data':[]}
            for file in request.FILES.getlist('documents'):
                save_file(file)
                file=str(file)
                print(file)
                if str(file).endswith('.sdf'):
                    file_path = MEDIA_ROOT + '/uploaded_data/'
                    with Chem.SDMolSupplier(file_path+file) as suppl:
                        for substance in suppl:
                            if substance is None:
                                continue
                            '''subst = Molecule()
                            subst.inchi_key = substance.GetProp('PUBCHEM_IUPAC_INCHIKEY')
                            subst.log_p = substance.GetProp('PUBCHEM_XLOGP3')
                            subst.num_h_acceptors = substance.GetProp('PUBCHEM_CACTVS_HBOND_ACCEPTOR')
                            subst.num_h_donors = substance.GetProp('PUBCHEM_CACTVS_HBOND_DONOR')
                            subst.molecular_mass = substance.GetProp('PUBCHEM_MOLECULAR_WEIGHT')
                            subst.save()'''
                    data['form']=form
                    data['existing_uploads']['data'].append([file,1,form['description'].value()])
        else:
            err_msg='Not a SDFile'
            data['err_msg']=err_msg
    return render(request,'upload.html',data)