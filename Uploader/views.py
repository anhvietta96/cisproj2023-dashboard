import os.path
from django.http import HttpRequest
from django.shortcuts import render
from .forms import SDFileForm, SDFileMult
from django.conf import settings
from compounds.file_handler import MoleculeIterator
from Uploader.file_handler import RequestFileIterator


def SDFile_Upload(request: HttpRequest):
    form = SDFileForm()
    data = {'form': form}
    if request.method == 'POST':
        form = SDFileForm(request.POST, request.FILES)
        if form.is_valid():
            form.save()
            file = form['document'].value()
            file_str = str(file)
            file_path = os.path.join(
                settings.MEDIA_ROOT, 'uploaded_data/', file_str)
            try:
                m = MoleculeIterator(file_path)
                m.iterate_over_molecules()
                data['form'] = form
                data['existing_uploads'] = {
                    'headers': ['File', 'Set', 'Description']}
                data['existing_uploads']['data'] = [
                    [file_str, 1, form['description'].value()]]
            except ValueError:
                data['err_msg'] = 'Not a valid SDFile'

        else:
            err_msg = 'Not a SDFile'
            data['err_msg'] = err_msg
    return render(request, 'Uploader/upload.html', data)


def SDFMultView(request: HttpRequest):
    form = SDFileMult()
    data = {'form': form}

    if request.method == 'POST':
        form = SDFileMult(request.POST, request.FILES)
        if form.is_valid():
            request_file_iterator = RequestFileIterator(
                os.path.join(settings.MEDIA_ROOT, 'uploaded_data/'),
                request.FILES.getlist('documents'),
                form)
            
            request_file_iterator.iterate_over_files()
            data = request_file_iterator.get_data_dict()
            request_file_iterator.add_to_set_from_form()
        else:
            data['form'] = form
            data['err_msg'] = 'An error occurred'
    return render(request, 'Uploader/upload.html', data)
