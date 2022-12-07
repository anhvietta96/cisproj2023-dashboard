import os.path
from django.shortcuts import render
from .forms import SDFileForm, SDFileMult
from dashboard.settings import MEDIA_ROOT
from compounds.file_handler import MoleculeIterator, FileIterator


class RequestFileIterator(FileIterator):
    def __init__(self, dirname, files):
        super().__init__(dirname)
        self.files = files

    def iterate_over_files(self):
        for file in self.files:
            file_path = os.path.join(self.dirname, file.name)

            # Save file
            with open(file_path, 'wb+') as destination:
                for chunk in file.chunks():
                    destination.write(chunk)

            super()._handle_mol_iterator(file_path)


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


def save_file(file, file_path):
    with open(file_path, 'wb+') as destination:
        for chunk in file.chunks():
            destination.write(chunk)


def SDFMultView(request):
    form = SDFileMult()
    data = {'form': form}

    if request.method == 'POST':
        form = SDFileMult(request.POST, request.FILES)
        if form.is_valid():
            data['form'] = form
            data['existing_uploads'] = \
                {'headers': ['File', 'Set', 'Description'], 'data': []}

            for file in request.FILES.getlist('documents'):
                file_path = os.path.join(
                    MEDIA_ROOT, 'uploaded_data/', file.name)
                save_file(file, file_path)
                try:
                    m = MoleculeIterator(file_path)
                    m.iterate_over_molecules()
                except ValueError:
                    data['err_msg'] = 'Not a valid SDFile'
                    continue

                data['existing_uploads']['data'].append(
                    [file, 1, form['description'].value()])
        else:
            err_msg = 'Not a valid SDFile'
            data['err_msg'] = err_msg
    return render(request, 'upload.html', data)
