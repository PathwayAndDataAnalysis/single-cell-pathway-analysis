from django.db import models
from mongoengine import Document, StringField


# Create your models here.

class UserData(Document):
    user_email = StringField()
    cell_data = StringField()

    def __str__(self):
        return self.user_email + ": " + self.cell_data
