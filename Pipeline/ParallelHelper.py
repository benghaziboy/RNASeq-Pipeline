from multiprocessing import Process, Pipe
from itertools import izip


def spawn(f):
    def fun(pipe, x):
        pipe.send(f(x))
        pipe.close()
    return fun


def parmap(f, X):
    pipe = [Pipe() for x in X]
    proc = [Process(target=spawn(f), args=(c, x))
            for x, (p, c) in izip(X, pipe)]
    [p.start() for p in proc]
    [p.join() for p in proc]
    return [p.recv() for (p, c) in pipe]


def pparmap(f, X, num_procs):
    res = []
    num_parts = int(len(X) / num_procs) + 1
    for i in range(num_parts):
        if ((i + 1) * num_procs) >= len(X):
            res += parmap(f, X[i * num_procs:])
        else:
            res += parmap(f, X[i * num_procs:((i + 1) * num_procs)])
    return res
