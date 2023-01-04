from django.db import models


# Create your models here.

class ThreadTask(models.Model):
    task_name = models.CharField(max_length=50, blank=True, null=True)
    is_done = models.BooleanField(blank=False, default=False)

    def __str__(self):
        return self.task_name
