from pyscipopt import Model, quicksum, SCIP_RESULT, SCIP_PARAMSETTING
from pyscipopt.linexpr import LinExpr
from data import read_data
from math import floor
from evalinstance import getAlpha
import timeit
import sys


def createModel(data):
    D, W, H, N, l, h, w = data['D'], data['W'], data['H'], data['N'], data['l'], data['h'], data['w']

    model = Model('compact2')
    V     = {}
    nodetype = {}
    E     = []

    # create all nodes
    s = 0
    nnodes = 1
    for i in reversed(N):
        V[i] = []
        for k in range(int(D / l[i])):
            V[i].append(nnodes)
            nodetype[nnodes] = i
            nnodes += 1

    # create all edges
    for (i,j) in [(i,j) for i in reversed(N) for j in reversed(N) if N.index(i) > N.index(j)]:
        pos = 0
        for v in V[i]:
            for k in range(int(l[i] / l[j])):
                E.append((v,V[j][pos]))
                pos += 1

    # connect every node with source
    for v in range(1, nnodes):
        E.append((s,v))

    x = {}
    y = {}
    p = {}

    for i in range(1, nnodes):
        p[i] = model.addVar('p_' + str(i), obj=0, lb=0, ub=H, vtype='C')

    for (i,j) in E:
        x[i,j] = model.addVar('x_' + str(i) + ',' + str(j), obj=-w[nodetype[j]], lb=0, ub=1, vtype='C')
        y[i,j] = model.addVar('y_' + str(i) + ',' + str(j), obj=0, lb=0, ub=1, vtype='B')

    # x <= y
    for (i,j) in E:
        model.addCons(x[i,j] - y[i,j] <= 0, name = 'bound_'+str(i)+'_'+str(j))

    # sum_{(i,j) in E : t(j) = n} x_{i,j} <= 1
    for n in N:
        model.addCons(quicksum(x[i,j] for (i,j) in E if nodetype[j] == n) <= 1, name='frac_' + str(n))

    # sum_{e in delta^{-}(i)} y_e <= 1
    for j in range(1, nnodes):
        model.addCons(quicksum(y[i,j] for (i,k) in E if k == j) <= 1, name='incoming_' + str(j))

    # sum_{e in delta^{+}(s) l[e] * y_e <= D}
    model.addCons(quicksum(l[nodetype[i]] * y[s,i] for i in range(1, nnodes)) <= D, name='depth')

    # sum_{e in delta^{+}(i)} l(e) * y_e <= sum_{e in delta^{-}(i)} l(e) * y_e
    for i in range(1, nnodes):

        # constraint is redundant for nodes without outgoing edges
        if nodetype[i] == N[0]:
            continue

        sum1 = quicksum(l[nodetype[j]] * y[i,j] for (k,j) in E if k == i)
        sum2 = quicksum(-l[nodetype[i]] * y[j,i] for (j,k) in E if k == i)
        model.addCons(sum1 + sum2 <= 0, name='length_' + str(i))

    # potential constraints with big-M
    for (i,j) in E:
        if i == 0:
            model.addCons(-H <= p[j] - h[nodetype[j]] * x[i,j] - H * y[i,j], name = 'potential_' + str(i) + '_' + str(j))
        else:
            model.addCons(-H <= p[j] - p[i] - h[nodetype[j]] * x[i,j] - H * y[i,j], name = 'potential_' + str(i) + '_' + str(j))

    # for j in range(1, nnodes):
    #     sum1 = quicksum(p[i]*y[i,k] for (i,k) in E if k == j if i != s)
    #     sum1 += quicksum(h[nodetype[k]] * x[i,k]  for (i,k) in E if k == j)
    #     sum1 += quicksum(-1 * p[k] for k in [j])
    #     cons = model.addCons(sum1 == 0, name='potential_' + str(j))

    return model, E, nodetype, x, y

def runCompact2(data, timelim, memlim, display, quite):
    """
    creates and solves instance with the reduced complete graph
    output:
    edgemap     - list containing all edge tuples (e,type,x,y)
    stats       - dict containing information about 'status', 'time', 'nnodes', 'readtime', 'nedges', 'primal', 'dual', 'dualroot', 'nxused', 'nyused'
    """
    model, E, nodetype, x, y = createModel(data)

    # set working limits
    model.setRealParam('limits/time', timelim)
    model.setRealParam('limits/memory', memlim)
    model.setHeuristics(SCIP_PARAMSETTING.AGGRESSIVE)

    # hide output
    if quite:
        model.hideOutput()

    model.optimize()

    if not quite:
        model.printBestSol()

    assert not model.getStatus() is 'infeasible' and not model.getStatus() is 'unbounded'

    # store statistics
    stats = {}
    stats['status']   = model.getStatus()
    stats['time']     = model.getSolvingTime()
    stats['nnodes']   = model.getNNodes()
    stats['dual']     = model.getDualbound()
    stats['dualroot'] = model.getDualboundRoot()
    stats['primal']   = model.getPrimalbound()
    stats['gap']      = model.getGap()
    stats['nxused']   = 0
    stats['nyused']   = 0
    stats['nedges']   = len(E)

    # compute edgemap
    edgemap = []
    nsols = len(model.getSols())
    assert nsols >= 0
    for e in E:
        t = nodetype[e[1]]
        if int(nsols) > 0:
            xval = float(model.getSolVal(model.getBestSol(), x[e]))
            yval = float(model.getSolVal(model.getBestSol(), y[e]))
            if yval >= 0.5:
                edgemap.append((e, t, xval, yval))
            if xval > 1e-06:
                stats['nxused'] += 1
            if yval > 1e-06:
                stats['nyused'] += 1
    edgemap.sort(key = lambda tup : tup[0][0])

    # display solution
    if display is True:
        from gui import displaySolution
        displaySolution(data, edgemap)

    return stats

# MAIN
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('usage: python compact2.py <instance> <time limit> <memory limit> [display] [quite]')
        exit(0)

    # read input
    data = read_data(str(sys.argv[1]))
    timelim = float(sys.argv[2])
    memlim = float(sys.argv[3])

    # optional parameter
    display = True if len(sys.argv) >= 5 and str(sys.argv[4]) == 'True' else False
    quite = True if len(sys.argv) >= 6 and str(sys.argv[5]) == 'True' else False

    # solve problem
    stats = runCompact2(data, timelim, memlim, display, quite)
    print('STATS: ' + str(stats))
