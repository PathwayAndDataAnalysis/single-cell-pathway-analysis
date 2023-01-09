from django.db import models
from mongoengine import StringField, Document, BooleanField


# Create your models here.


class Analysis(Document):
    analysisName = StringField(required=True)
    isFilteringDone = models.BooleanField(default=False)
    isPCADone = models.BooleanField(default=False)
    isUMAPDone = models.BooleanField(default=False)
    isDone = BooleanField(default=False)

    def __str__(self):
        return self.analysisName
