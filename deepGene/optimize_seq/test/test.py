def throw():
    try:
        raise IndexError('index error in throw')
    except Exception:
        pass


def my_function():
    try:
        throw()
    except Exception, e:
        print e


if __name__ == '__main__':
    my_function()
