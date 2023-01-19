from django.http import HttpResponse


# Create your views here.

def httpMethodError(expectedMethod: str, receivedMethod: str):
    """
    Returns an HttpResponse with a 404 status code and a message indicating that the expected method was expectedMethod
    @param expectedMethod:
    @param receivedMethod: str
    @return: HttpResponse
    """
    print('Expected method: ' + expectedMethod + ', but received method: ' + receivedMethod)
    return HttpResponse({'status': 'ERROR',
                         'message': 'Expected method: ' + expectedMethod + ', but received method: ' + receivedMethod},
                        status=404)


def exceptionInRequest(e: Exception):
    """
    Returns an HttpResponse with a 400 status code and a message indicating that an exception was raised
    @param e: Exception
    @return: HttpResponse
    """
    print("Exception in request: " + str(e))
    return HttpResponse({'status': 'ERROR',
                         'message': 'Exception occurred: ' + str(e)}, status=400)
