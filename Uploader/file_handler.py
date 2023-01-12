import os.path
from compounds.file_handler import MoleculeIterator
from django.forms import Form


class RequestFileIterator(MoleculeIterator):
    """
    Class for the iteration over all files uploaded by the user
    """

    def __init__(self, dirname: str, files, form: Form):
        """
        :param dirname: dir_path to the directory of the uploaded files
        :param files: uploaded files from the user
        :type files: HttpRequest.FILES
        :param form: validated upload form
        """
        if not os.path.exists(dirname):
            os.makedirs(dirname)

        super().__init__()
        self.dirname = dirname
        self.files = files
        self.form = form
        self.set_name = self.form['set_name'].value()
        self.set_description = self.form['description'].value()

        # For Front-End communication
        self.data = {
            'form': self.form,
            'existing_uploads': {
                'headers': ['File', 'Set', 'Description'],
                'data': []}
        }

    def iterate_over_files(self) -> None:
        """
        Function for the iteration over all uploaded files by the user.
        Molecules will be saved in the Molecule model
        """
        for file in self.files:
            file_path = os.path.join(self.dirname, file.name)

            # Save file
            with open(file_path, 'wb+') as destination:
                for chunk in file.chunks():
                    destination.write(chunk)
            if super().iterate_over_molecules(file_path):
                self.data['existing_uploads']['data'].append(
                    [file, self.set_name, self.set_description])

    def get_data_dict(self) -> dict:
        """
        Return the data dictionary for front-end communication
        :return: data dictionary containing information about the form,
            uploaded files and errors
        """
        self.data['err_msg'] = "\n".join(self.err_msgs)
        return self.data

    def add_to_set_from_form(self) -> None:
        """
        Adds the molecules to the set from the given form
        """
        super().add_to_set(self.set_name, self.set_description)
