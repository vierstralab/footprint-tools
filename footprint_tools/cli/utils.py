from argh.exceptions import CommandError

def tuple_ints(arg):
    """
    Function to parser a tuple of integers from command line
    """
    try:
        fw, rev = list(map(int, arg.split(',')))
        return (fw, rev)
    except Exception as e:
        raise CommandError("Argument must be a comma-delimited tuple -- i.e., 0,-1")

def list_ints(arg):
    """
    Function to parser a list of integers from command line
    """
    try:
        return list(map(int, arg.split(',')))
    except Exception as e:
        raise CommandError("Argument must be a comma-delimited list -- i.e., 0,1,2,3")

def get_kwargs(keys, kwargs):
    return {k:kwargs[k] for k in keys if k in kwargs}

def chunkify(l, n):
    """Splits an iterable list into n chunks"""
    return [l[i::n] for i in range(n)]
