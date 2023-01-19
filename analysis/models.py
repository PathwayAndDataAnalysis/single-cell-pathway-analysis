import enum

from mongoengine import StringField, Document, IntField


# Create your models here.

class TASK_STATUS(enum.Enum):
    DEFAULT = 0
    COMPLETED = 1
    FAILED = -1
    RUNNING = 2


class Analysis(Document):
    analysisName = StringField(required=True)
    isFilteringDone = IntField(default=TASK_STATUS.DEFAULT.value, required=True, min_value=-1, max_value=2)
    isPCADone = IntField(default=TASK_STATUS.DEFAULT.value, min_value=-1, max_value=2)
    isUMAPDone = IntField(default=TASK_STATUS.DEFAULT.value, min_value=-1, max_value=2)
    isAllDone = IntField(default=TASK_STATUS.DEFAULT.value, min_value=-1, max_value=2)
    errorMessage = StringField(default="")

    def __str__(self):
        return self.analysisName
