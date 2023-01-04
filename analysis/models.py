from django.db import models
from mongoengine import StringField, Document, BooleanField


# Create your models here.


class Analysis(Document):
    analysis_name = StringField()
    # is_filtering_done = models.BooleanField(default=False)
    # is_pca_done = models.BooleanField(default=False)
    # is_umap_done = models.BooleanField(default=False)
    isDone = BooleanField(default=False)

    def __str__(self):
        return self.analysis_name
