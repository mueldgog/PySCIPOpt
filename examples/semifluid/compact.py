from pyscipopt import Model, SCIP_RESULT, SCIP_PARAMSETTING
from data import read_data
from math import floor
from evalinstance import getAlpha
import timeit
import sys

def createModel(data, readtimelim):
    """creates and returns the model which uses the full graph"""

    # input
    D, W, H, N, l, h, w = data['D'], data['W'], data['H'], data['N'], data['l'], data['h'], data['w']

    model = Model('compact')

    E = []
    F = [0]
    T = {}
    nnodes = 1
    cons_length = {}
    cons_frac   = {}
    xvar = {}
    yvar = {}
    pvar = {}

    # compute dominance relations
    alpha = getAlpha(D,l,N)

    # create compact model
    start = timeit.default_timer()
    while len(F) > 0 and timeit.default_timer() - start < readtimelim:
        u = F.pop()

        lu = D if u == 0 else l[T[u]]
        items = N if u == 0 else [item for item in N if N.index(item) < N.index(T[u])]

        for t in items:
            lt = l[t]

            n = alpha[t,-1] if u == 0 else alpha[t,T[u]]

            for k in range(n):
                # add edge to the graph
                e = (u,nnodes)
                nnodes += 1

                E.append((e,t))

                F.append(e[1])
                T[e[1]] = t

                # create variables
                x = model.addVar('x_'+str(e[0])+'_'+str(e[1]), obj=-w[t], lb=0, ub=1, vtype='C')
                y = model.addVar('y_'+str(e[0])+'_'+str(e[1]), obj=0, lb=0, ub=1, vtype='B')
                p = model.addVar('p_'+str(e[1]), obj=0, lb=0, ub=H, vtype='C')
                xvar[e] = x
                yvar[e] = y
                pvar[e[1]] = p

                # create constraints
                model.addCons(x - y <= 0, name = 'bound_'+str(e[0])+'_'+str(e[1]))

                if u == 0:
                    model.addCons(p - h[t]*x == 0, name = 'potential_'+str(e[0])+'_'+str(e[1]))
                else:
                    model.addCons(p - pvar[u] - h[t]*x == 0, name = 'potential_'+str(e[0])+'_'+str(e[1]))

                if not t in cons_frac:
                    cons_frac[t] = model.addCons(x <= 1, name = 'frac_'+str(t))
                else:
                    model.addConsCoeff(cons_frac[t], x, 1)

                if not 0 in cons_length:
                    cons_length[0] = model.addCons(lt*y <= D, name = 'length_'+str(e[0]))
                else:
                    model.addConsCoeff(cons_length[e[0]], y, lt)

                cons_length[e[1]] = model.addCons(-lt*y <= 0, name = 'length_'+str(e[1]))

    stop = timeit.default_timer() - start
    return model, stop, E, xvar, yvar


def runCompact(data, timelim, memlim, readtimelim, solfileread, solfilewrite, display, quite):
    """
    creates and solves instance with the complete graph
    output:
    edgemap     - list containing all edge tuples (e,type,x,y)
    stats       - dict containing information about 'status', 'time', 'nnodes', 'readtime', 'nedges', 'primal', 'dual', 'dualroot', 'nxused', 'nyused'
    """
    stats = {}
    edgemap = []

    # initialize status
    stats['status'] = 'unknown'
    stats['time'] =  0.0
    stats['nnodes'] = 0
    stats['readtime'] = 0.0
    stats['nedges'] = 0
    stats['nxused'] = 0
    stats['nyused'] = 0

    # note that E is a list of tuples (e,idx,type)
    model, readtime, E, xvars, yvars = createModel(data, readtimelim)
    assert (not model is None) and (not E is None) and not xvars is None and (not yvars is None)

    stats['readtime'] = readtime
    stats['primal']   = model.infinity()
    stats['dual']     = -model.infinity()
    stats['dualroot'] = -model.infinity()
    stats['gap']      = model.infinity()
    stats['nedges']   = len(E)

    # stop if we reached the time limit
    if float(readtime) >= float(readtimelim):
        stats['status'] = 'readtime'
        return edgemap, stats

    # hide output
    if quite:
        model.hideOutput()

    # read solution file
    if not solfileread is None:
        model.readSol(str(solfileread))

    # set working limits
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.setSeparating(SCIP_PARAMSETTING.OFF)
    model.setHeuristics(SCIP_PARAMSETTING.AGGRESSIVE)
    model.setRealParam('limits/time', timelim)
    model.setRealParam('limits/memory', memlim)

    model.optimize()
    assert not model.getStatus() is 'infeasible' and not model.getStatus() is 'unbounded'

    # update statistics
    stats['status']   = model.getStatus()
    stats['time']     = model.getSolvingTime()
    stats['nnodes']   = model.getNNodes()
    stats['dual']     = model.getDualbound()
    stats['dualroot'] = model.getDualboundRoot()
    stats['primal']   = model.getPrimalbound()
    stats['gap']      = model.getGap()
    stats['nxused']   = 0
    stats['nyused']   = 0

    # compute edgemap
    nsols = len(model.getSols())
    assert nsols >= 0
    for (e,t) in E:
        if int(nsols) == 0:
            edgemap.append((e, t, 0.0, 0.0))
        else:
            xval = float(model.getVal(xvars[e]))
            yval = float(model.getVal(yvars[e]))
            edgemap.append((e, t, xval, yval))
            if xval > 1e-06:
                stats['nxused'] += 1
            if yval > 1e-06:
                stats['nyused'] += 1

    # write solution file
    if not solfilewrite is None:
        model.writeBestSol(str(solfilewrite))

    # display solution
    if display is True:
        from gui import displaySolution
        displaySolution(data, edgemap)

    return edgemap, stats


# MAIN
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('usage: python compact.py <instance> <time limit> <memory limit> [display] [quite]')
        exit(0)

    # read input
    data = read_data(str(sys.argv[1]))
    timelim = float(sys.argv[2])
    memlim = float(sys.argv[3])

    # optional parameter
    display = True if len(sys.argv) >= 5 and str(sys.argv[4]) == 'True' else False
    quite = True if len(sys.argv) >= 6 and str(sys.argv[5]) == 'True' else False

    edgemap, stats = runCompact(data, timelim, memlim, timelim, None, None, display, quite)

    # print statistics
    print('STATS: ' + str(stats))
