import random, gzip
import sys

if sys.version_info >= (3, 0):
   from queue import Queue
else:
   from Queue import Queue

from fractions import Fraction
from decimal import Decimal

def read_data(filename):
    """read_data(filename): reads data file
    returns:
        * D,W,H: container's width, height, depth
        * n: number of semifluids to put in the data
        * l: dict of lengths of items
        * v:  ... volumes (fraction of total volume)
        * w:  ... unit values
        """
    from collections import deque
    data = deque()
    # with gzip.open(filename, "rt", encoding="ascii") as f:
    with open(filename, "rt") as f:
        for line in f:
            line += "#"
            nocomment = line[:line.find("#")]
            data += nocomment.split()
    D = Fraction(data.popleft())
    W = Fraction(data.popleft())
    H = Fraction(data.popleft())
    n = int(data.popleft())

    l = {}    # lengths
    v = {}    # volumes (fraction of total volume)
    w = {}    # unit values
    h = {}
    for i in range(n):
        l[i] = Fraction(data.popleft())
        v[i] = Fraction(data.popleft())
        w[i] = Fraction(data.popleft()) * v[i]
        h[i] = v[i] * H * D / l[i]

    # sort semifluids w.r.t. lengths
    N = sorted(l, key=l.__getitem__)

    return {'D': D, 'W':W, 'H':H, 'N':N, 'l':l, 'h':h, 'w':w, 'v':v}


# def read_data(filename):
#    """Reads data for semi-fluid packing instance from a given file.

#    Keyword arguments:
#    filename --- the name of the file to be used
#    Returns:
#       * D,W,H: container's depth, width, height
#       * N: list of all semifluids sorted w.r.t. their length
#       * l: dict of lengths of items
#       * h: dict of (maximum) heights of items
#       * w: dict of unit values
#    """
#    f = open(filename, 'r')

#    D = -1.0
#    W = -1.0
#    H = -1.0
#    l = {}
#    h = {}
#    w = {}

#    if sys.version_info >= (3, 0):
#       q = Queue()
#    else:
#       q = Queue()


#    for line in f:
#       # ignore comments
#       if line[0] == '#':
#          continue
#       q.put(line)
#    f.close()

#    D = Fraction(q.get())
#    W = Fraction(q.get())
#    H = Fraction(q.get())
#    n = Fraction(q.get())

#    for i in range(int(n)):
#        split = q.get().split()
#        l[i] = Fraction(split[0])
#        h[i] = Fraction(split[1]) * H * D / l[i]
#        w[i] = Fraction(split[2])

#    # sort semifluids w.r.t. lengths
#    N = sorted(l, key=l.__getitem__)

#    return {'D': D, 'W':W, 'H':H, 'N':N, 'l':l, 'h':h, 'w':w}
#    D,W,H,n,l,v,w
