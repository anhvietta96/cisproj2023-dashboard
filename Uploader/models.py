from django.db import models
from django import forms
# Create your models here.

class SDFile(models.Model):
    description = models.CharField(max_length=100,blank=True)
    document=models.FileField(upload_to='uploaded_data/')
    uploaded_at=models.DateTimeField(auto_now_add=True)

    def path():
        return '/media/uploaded_data'