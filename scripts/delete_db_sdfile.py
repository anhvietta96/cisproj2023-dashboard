from Uploader.models import SDFile

def run():
    for file in SDFile.objects.all().iterator():
        file.delete()