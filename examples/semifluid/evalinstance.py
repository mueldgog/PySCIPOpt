#########################################################################
#                                                                       #
#                  This file is part of the program and library         #
#         SCIP --- Solving Constraint Integer Programs                  #
#                                                                       #
#    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                        #
#                          fuer Informationstechnik Berlin              #
#                                                                       #
#  SCIP is distributed under the terms of the ZIB Academic License.     #
#                                                                       #
#  You should have received a copy of the ZIB Academic License.         #
#  along with SCIP; see the file COPYING. If not email to scip@zib.de.  #
#                                                                       #
#########################################################################

from data import read_data
from math import floor
from os import listdir
import sys

def getAlpha(D,l,N):
    # alpha[i,k] = min{l[j] / l[i] - 1} for all i < j < k with l[j] % l[i] = 0
    alpha = {}
    for (i,k) in [(i,k) for i in N for k in N if N.index(i) < N.index(k)]:
        alpha[i,k] = int(l[k] / l[i])
        for j in [j for j in N if N.index(i) < N.index(j) and N.index(j) < N.index(k) and (l[j] % l[i]) == 0]:
            alpha[i,k] = min(alpha[i,k], int(l[j] / l[i] - 1))
        assert alpha[i,k] >= 0

    if len(N) >= 3:
        i,j,k = N[0], N[1], N[2]
        if int(l[k] / l[i]) == int(l[k] / l[j]) * int(l[j] / l[i]) + int((l[k] % l[j]) / l[i]):
            alpha[i,k] = min(alpha[i,k], int((l[k] % l[j]) / l[i]))

    for i in N:
        alpha[i,-1] = int(D / l[i])
        for j in [j for j in N if N.index(i) < N.index(j) and (l[j] % l[i]) == 0]:
            alpha[i,-1] = min(alpha[i,-1], int(l[j] / l[i] - 1))
        assert alpha[i,-1] >= 0

    return alpha

def computeNEdgesComplex(D, l, N):
    """Recursive formular to compute the total number of edges when using symmetry breaking."""
    alpha = getAlpha(D,l,N)
    T = {N[0] : 1}
    for pos1 in range(1,len(N)):
        val = 0
        t1 = N[pos1]
        for pos2 in range(0,pos1):
            t2 = N[pos2]
            val += T[t2]*alpha[t2,t1]
        T[N[pos1]] = 1 + val

    res = 0
    for t in N:
        res += T[t] * alpha[t,-1]

    return res


def computeNEdgesTrivial(D, l, N):
    """Recursive formular to compute the total number of edges."""
    alpha = {}
    for (t1,t2) in [(x,y) for x in N for y in N if N.index(x) < N.index(y)]:
        alpha[t1,t2] = int(floor(l[t2] / l[t1]))
    for t in N:
        alpha[t,-1] = int(floor(D / l[t]))

    T = {N[0] : 1}
    for pos1 in range(1,len(N)):
        val = 0
        t1 = N[pos1]
        for pos2 in range(0,pos1):
            t2 = N[pos2]
            val += T[t2]*alpha[t2,t1]
        T[N[pos1]] = 1 + val

    res = 0
    for t in N:
        res += T[t] * alpha[t,-1]

    return res


def computeNEdgesComplex2(data):
    D, W, H, N, l, h, w = data['D'], data['W'], data['H'], data['N'], data['l'], data['h'], data['w']

#    model = Model('compact2')
    V     = [-1]
    vtype = {}
    E     = []

    # create all nodes
    for i in N:
        for k in range(int(D / l[i])):
            v = len(V) - 1
            V.append(v)
            vtype[v] = i

    # create all edges
    nedges = 0
    alpha = getAlpha(D,l,N)

    for i in range(len(N) + 1):
        t1 = -1 if i == len(N) else N[i]
        l1 = D  if t1 == -1 else l[t1]
        for j in range(i):
            t2 = N[j]
            nedges += int(l1 / l[t2]) * int(D / l1)

    return nedges


if __name__ == '__main__':

    if len(sys.argv) != 2:
        instances = listdir('./data/')
        instances.sort()

        for instance in instances:
            if not 'hard' in instance:
                continue
            data = read_data('./data/' + str(instance))
            val1 = computeNEdgesComplex2(data)
            val2 = computeNEdgesComplex(data['D'], data['l'], data['N'])
#            assert val1 >= val2
            print('%20s %10e %10e   FACTOR: %10f' % (instance, float(val1), float(val2), float(val2) / float(val1)))
    else:
        data = read_data(str(sys.argv[1]))
        val1 = computeNEdgesComplex2(data)
        val2 = computeNEdgesComplex(data['D'], data['l'], data['N'])
#        assert val1 >= val2
        print(str(sys.argv[1]), val1, val2)

