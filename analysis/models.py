from mongoengine import StringField, Document, BooleanField

# Create your models here.

class Analysis(Document):
    analysisName = StringField(required=True)
    isFilteringDone = BooleanField(default=False)
    isPCADone = BooleanField(default=False)
    isUMAPDone = BooleanField(default=False)
    isAllDone = BooleanField(default=False)

    def __str__(self):
        return self.analysisName
