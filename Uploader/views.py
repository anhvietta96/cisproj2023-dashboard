import os.path
from django.shortcuts import render
from .forms import SDFileForm, SDFileMult
from dashboard.settings import MEDIA_ROOT
from compounds.file_handler import MoleculeIterator
from Uploader.file_handler import RequestFileIterator


def SDFile_Upload(request):
    form = SDFileForm()
    data = {'form': form}
    if request.method == 'POST':
        form = SDFileForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            file = form['document'].value()
            file = str(file)
            file_path = os.path.join(MEDIA_ROOT, 'uploaded_data/', file)
            try:
                m = MoleculeIterator(file_path)
                m.iterate_over_molecules()
                data['form'] = form
                data['existing_uploads'] = {
                    'headers': ['File', 'Set', 'Description']}
                data['existing_uploads']['data'] = [
                    [file, 1, form['description'].value()]]
            except ValueError:
                data['err_msg'] = 'Not a valid SDFile'

        else:
            err_msg = 'Not a SDFile'
            data['err_msg'] = err_msg
    return render(request, 'upload.html', data)


def SDFMultView(request):
    form = SDFileMult()
    data = {'form': form}

    if request.method == 'POST':
        form = SDFileMult(request.POST, request.FILES)
        if form.is_valid():
            request_file_iterator = RequestFileIterator(
                os.path.join(MEDIA_ROOT, 'uploaded_data/'),
                request.FILES.getlist('documents'),
                form)
            request_file_iterator.iterate_over_files()
            data = request_file_iterator.get_data_dict()
            request_file_iterator.add_to_new_set()
        else:
            data['err_msg'] = 'Not a valid SDFile'
    return render(request, 'upload.html', data)
