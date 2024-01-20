import warnings
import inspect


def chunks(lst, n):
    """
    Yield successive n-sized chunks from lst.
    https://stackoverflow.com/a/312464/17572887
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def verbose(f):
    """
    Decorator to run functions in silenced/verbose modes.
    """
    def wrapper(*args, **kwargs):
        flag = True
        if 'verbose' in kwargs:
            flag = kwargs['verbose']
            if 'verbose' not in inspect.signature(f).parameters:
                del kwargs['verbose']

        if flag:
            return f(*args, **kwargs)
        else:
            with warnings.catch_warnings():
                # Suppress warnings
                warnings.simplefilter('ignore')
                return f(*args, **kwargs)

    return wrapper
