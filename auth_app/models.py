# Create your models here.
from mongoengine import *


class Users(Document):
    name = StringField()
    email = StringField()
    password = StringField()

    def __str__(self):
        return self.name + " " + self.email + " " + self.password
