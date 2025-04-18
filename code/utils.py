import yaml
from collections import namedtuple
from time import time

    
def timer(func):
    def wrapper(*args,**kw):
        start = time()
        res = func(*args,**kw)
        duration = time() - start
        print(f"run {func.__name__} in {duration:.1f} seconds")
        return res
    return wrapper
    
def yaml2dict(path):
    with open(path) as f:
        x = yaml.safe_load(f)
    res = {}
    for i in x:
        res[i] = x[i]['value']
    return res

def dict2namedtuple(dic):
    return namedtuple('Config', dic.keys())(**dic)

def load_yaml(path):
    res = yaml2dict(path)
    config = dict2namedtuple(res)
    print(config)
    return config

if __name__ == '__main__':
    with open('p2test.yaml') as f:
        x = yaml.safe_load(f)
    print(x)
    res = {}
    for i in x:
        res[i] = x[i]['value']
        print(res[i])
    config = dict2namedtuple(res)
    print(config)
