import os.path
from compounds.file_handler import FileIterator
from compounds.models import MoleculeSet
from time import perf_counter


class RequestFileIterator(FileIterator):
    def __init__(self, dirname, files, form):
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        super().__init__(dirname)
        self.files = files
        self.form = form
        self.set_name = self.form['set_name'].value()
        self.set_description = self.form['description'].value()

        # Front-End
        self.data = {
            'form': self.form,
            'existing_uploads': {
                'headers': ['File', 'Set', 'Description'],
                'data': []}
        }

    def iterate_over_files(self):
        for file in self.files:
            file_path = os.path.join(self.dirname, file.name)

            # Save file
            with open(file_path, 'wb+') as destination:
                for chunk in file.chunks():
                    destination.write(chunk)

            if super()._handle_mol_iterator(file_path):
                self.data['existing_uploads']['data'].append(
                    [file, self.set_name, self.set_description])

    def get_data_dict(self):
        self.data['err_msg'] = "\n".join(self.err_msgs)
        return self.data

    def add_to_new_set(self):
        mol_set = MoleculeSet(
            set_name=self.set_name,
            description=self.set_description)
        mol_set.save()
        mol_set.molecules.add(*self.get_mol_list())
        mol_set.save()
