from django.db import models
from django.core.validators import FileExtensionValidator
from dashboard.settings import MEDIA_ROOT
import os

class SDFile(models.Model):
    description = models.CharField(max_length=100, blank=True)
    document = models.FileField(
        upload_to='uploaded_data/',
        validators=[FileExtensionValidator(allowed_extensions=['sdf'])])
    uploaded_at = models.DateTimeField(auto_now_add=True)

    def path(self):
        return os.path.join(MEDIA_ROOT,'uploaded_data')
