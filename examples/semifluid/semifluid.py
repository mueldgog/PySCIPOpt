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

from pyscipopt import Model, Pricer, LP, SCIP_RESULT, SCIP_PARAMSETTING
from data import read_data
from evalinstance import getAlpha
from math import floor
from gui import displaySolution
import sys


def DEBUGMSG(msg):
    return
    print(msg)


class fluidPricer(Pricer):

    def updateCols(self):
        """transfers new bounds of SCIP variables to their corresponding column bounds"""
        lp = self.data['LP']
        for var in list(self.data['var2col']):
            col = self.data['var2col'][var]
            lp.chgBound(col, var.getLbLocal(), var.getUbLocal())


    def pricerredcost(self):
        """the reduced cost function for the variable pricer"""

        # stop if there is no front edge remaining
        if len(self.data['F']) == 0:
            return {'result':SCIP_RESULT.SUCCESS}

        # update columns
        self.updateCols()

        lp = self.data['LP']
        l  = self.data['data']['l']
        N  = self.data['data']['N']

        # fix frontier variables
        for (u,v) in self.data['F']:
            x_col = self.data['LP_edge2x'][(u,v)]
            y_col = self.data['LP_edge2y'][(u,v)]
            lp.chgBound(x_col, 0, 0)
#            lp.chgBound(y_col, 0, 0)

        # first solve the problem with fixed frontier variables
        lp.solve()

        # unfix frontier variables
        for (u,v) in self.data['F']:
            x_col = self.data['LP_edge2x'][(u,v)]
            y_col = self.data['LP_edge2y'][(u,v)]
            lp.chgBound(x_col, 0, 1)
#            lp.chgBound(y_col, 0, 1)

        # solve LP with iteration limit; collect basis indicies before and after the LP solve
        basis_pre = lp.getBasisInds()
        lp.solve(True)#, int(self.data['lpiter']))
        imprcols = list(set(lp.getBasisInds()) - set(basis_pre))
        redcosts = lp.getRedcost()


#        print("----------------------")
#        primal = lp.getPrimal()

        # update frontier
        addF = []
        for (u,v) in self.data['F']:
            type_uv = self.data['LP_edge2item'][(u,v)]
            x_col   = self.data['LP_edge2x'][(u,v)]
            y_col   = self.data['LP_edge2y'][(u,v)]
#            print((u,v),redcosts[x_col], redcosts[y_col], primal[x_col], primal[y_col])

            if redcosts[x_col] + redcosts[y_col] <= -1e-08 or x_col in imprcols or y_col in imprcols:
                DEBUGMSG('%s redcost = %e' % (str(x_col), redcosts[x_col]))
                DEBUGMSG('%s in basis = %d' % (str(x_col), x_col in imprcols))
                DEBUGMSG('%s redcost = %e' % (str(y_col), redcosts[y_col]))
                DEBUGMSG('%s in basis = %d' % (str(y_col), y_col in imprcols))
                DEBUGMSG('remove edge (%d, %d)' % (u, v))

                self.data['F'].remove((u,v))
                addEdgeSCIP(self.model,self.data, u, v)
#                print("addEdgeSCIP ", u, v, type_uv)
                lv = l[type_uv]

                # update the frontier; consider all items which fit on top of item type[u,v]
                for idx in range(N.index(type_uv)):
                    t = N[idx]
                    lt = self.data['data']['l'][t]

                    for k in range(self.data['alpha'][t,type_uv]):
                        DEBUGMSG('add edge: ' + str((v, self.data['LP_nnodes'])) + ' type =' + str(t))
                        e = addEdgeLP(self.data, v, t)
#                        print("addEdgeLP ", e, t)
                        addF.append(e)


        # add new edges to the frontier
        self.data['F'] += addF

        return {'result':SCIP_RESULT.SUCCESS}


    # The initialisation function for the variable pricer to retrieve the transformed constraints of the problem
    def pricerinit(self):

        # transform all existing constraints
        for e in list(self.data['cons_length']):
            self.data['cons_length'][e] = self.model.getTransformedCons(self.data['cons_length'][e])
        for e in list(self.data['cons_frac']):
            self.data['cons_frac'][e] = self.model.getTransformedCons(self.data['cons_frac'][e])
        return

def addEdgeLP(pricerdata, u, t):
    """Adds an edge to the LP of the pricer. The corresponding x and y variables are fixed to zero.

    Keyword arguments:
    data --- pricer data
    u    --- starting node of the edge
    t    --- item type of the edge

    Returns:
    tuple (u,v) where v is the new created node
    """
    lp = pricerdata['LP']
    v  = int(pricerdata['LP_nnodes'])
    ht = float(pricerdata['data']['h'][t])
    wt = float(pricerdata['data']['w'][t])
    lt = float(pricerdata['data']['l'][t])
    edge2lengthrow = pricerdata['LP_edge2lengthrow']

    # add x_{u,v} to (frac) constraint
    xfracrow = pricerdata['LP_item2fracrow'][t]

    # add l[u,v] * y_{u,v] to predecessor (length) constraints
    ylengthrow = edge2lengthrow[pricerdata['LP_node2edge'][u]]

    pv_col = lp.ncols()
    lp.addCol([], 0, 0, pricerdata['data']['H'])
    x_col = lp.ncols()
    lp.addCol([(xfracrow,1)], -wt, 0, 0)
    y_col = lp.ncols()
    lp.addCol([(ylengthrow,lt)], 0, 0, 0)

    entries = []
    lhss = []
    rhss = []

    # add (length) constraint for edge (u,v)
    edge2lengthrow[u,v] = lp.nrows()
    entries.append([(y_col, -lt)])
    lhss.append(-lp.infinity())
    rhss.append(0)

    # add p_v - p_u - x_{u,v} * h_t = 0
    if u == 0:
        entries.append([(pv_col, 1), (x_col, -ht)])
        lhss.append(0)
        rhss.append(0)
    else:
        pu_col = pricerdata['LP_node2p'][u]
        entries.append([(pv_col, 1), (pu_col, -1), (x_col, -ht)])
        lhss.append(0)
        rhss.append(0)

    # add x_{u,v} - y_{u,v} <= 0
    entries.append([(x_col, 1), (y_col, -1)])
    lhss.append(-lp.infinity())
    rhss.append(0)

    # add all rows at once
    lp.addRows(entries, lhss, rhss)

    # update pricer data structure
    pricerdata['LP_edge2x'][(u,v)]    = x_col
    pricerdata['LP_edge2y'][(u,v)]    = y_col
    pricerdata['LP_node2p'][v]        = pv_col
    pricerdata['LP_nnodes']           += 1
    pricerdata['LP_edge2item'][(u,v)] = t

    # set predecessor of edge of node v
    pricerdata['LP_node2edge'][v] = (u,v)
    return (u,v)


# adds all variables and coefficients to SCIP
def addEdgeSCIP(model, pricerdata, u, v, pricedvar=True):
    D = pricerdata['data']['D']
    t = pricerdata['LP_edge2item'][(u,v)]
    lt = pricerdata['data']['l'][t]
    ht = pricerdata['data']['h'][t]
    xobj = -float(pricerdata['data']['w'][t])

    # create p_v, x_{u,v}, and y_{u,v}
    p  = model.addVar('p_'+str(v), obj=0, lb=0, ub=pricerdata['data']['H'], vtype='C', pricedVar=pricedvar)
    x  = model.addVar('x_'+str(u)+'_'+str(v), obj=xobj, lb=0, ub=1, vtype='C', pricedVar=pricedvar)
    y  = model.addVar('y_'+str(u)+'_'+str(v), obj=0, lb=0, ub=1, vtype='B', pricedVar=pricedvar)

    pricerdata['var_p'][v] = p
    pricerdata['var_x'][(u,v)] = x
    pricerdata['var_y'][(u,v)] = y

    # store mapping of each variable to it's corresponding LP column
    pricerdata['var2col'][p] = pricerdata['LP_node2p'][v]
    pricerdata['var2col'][x] = pricerdata['LP_edge2x'][(u,v)]
    pricerdata['var2col'][y] = pricerdata['LP_edge2y'][(u,v)]

    # add x_u,v <= y_u,v
    model.addCons(x - y <= 0, name = 'bound_'+str(u)+'_'+str(v))

    # create corresponding length constraint or add coefficient
    if not 0 in pricerdata['cons_length']:
        pricerdata['cons_length'][0] = model.addCons(lt*y <= D, name = 'length_'+str(u), separate = False, modifiable = True)
    else:
        assert u in pricerdata['cons_length']
        model.addConsCoeff(pricerdata['cons_length'][u], y, lt)
    pricerdata['cons_length'][v] = model.addCons(-lt*y <= 0, name = 'length_'+str(v), separate = False, modifiable = True)

    # add p_v = p_u  + x_{u,v}*h_{u,v}
    if u == 0:
        model.addCons(1 * p - ht * x == 0, name = 'potential_'+str(u)+'_'+str(v))
    else:
        model.addCons(1 * p - 1 * pricerdata['var_p'][u] - ht * x == 0, name = 'potential_'+str(u)+'_'+str(v))

    # add coefficient to (frac) constraint
    if not t in pricerdata['cons_frac']:
        pricerdata['cons_frac'][t] = model.addCons(x <= 1, name = 'frac_'+str(t), separate = False, modifiable = True)
    else:
        model.addConsCoeff(pricerdata['cons_frac'][t], x, 1)


def createModel(model, pricer, data, soledgemap, lpiter):
    """creates and returns the model which uses some starting edges"""
    D = data['D']
    l = data['l']
    h = data['h']
    N = data['N']
    alpha = getAlpha(D,l,N)

    # creating a pricer
    pricer.data = {}
    pricer.data['data']       = data
    pricer.data['soledgemap'] = soledgemap
    pricer.data['F']          = []
    pricer.data['lpiter']     = lpiter
    pricer.data['alpha']      = alpha

    # LP data
    pricer.data['LP']                = LP('LP pricer')
    pricer.data['LP_item2fracrow']   = {} # maps each item to it's corresponding (frac) constraint
    pricer.data['LP_edge2lengthrow'] = {} # maps each edge to it's corresponding (length) constraint
    pricer.data['LP_edge2x']         = {} # maps each edge to it's corresponding x column
    pricer.data['LP_edge2y']         = {} # maps each edge to it's corresponding y column
    pricer.data['LP_edge2item']      = {} # maps each edge to it's corresponding item type
    pricer.data['LP_node2edge']      = {0 : (-1,0)}
    pricer.data['LP_node2p']         = {}
    pricer.data['LP_nnodes']         = 1

    # create empty (frac) LP rows
    for i in N:
        pricer.data['LP_item2fracrow'][i] = i
        pricer.data['LP'].addRow([], -pricer.data['LP'].infinity(), 1)

    # create empty (length) LP rows
    pricer.data['LP_edge2lengthrow'] = {(-1,0) : pricer.data['LP'].nrows()}
    pricer.data['LP'].addRow([], -pricer.data['LP'].infinity(), D)

    # SCIP data
    pricer.data['var_p'] = {}
    pricer.data['var_x'] = {}
    pricer.data['var_y'] = {}
    pricer.data['var2col'] = {}
    pricer.data['cons_length'] = {}
    pricer.data['cons_frac']   = {}

    # add starting frontier; depends on the start edges
    if soledgemap is None:
        for item in N:
            for k in range(pricer.data['alpha'][item,-1]):
                e = addEdgeLP(pricer.data, 0, item)
#                print("addEdgeLP ", e, item)
                pricer.data['F'].append(e)

    else:
        itemOnNode = {}
        node2type = {0 : -1}
        nodemap = {}
        idx = 0
        maxp = {0 : 0}

        for (e,item,xval,yval) in soledgemap:
            if float(yval) < 0.5:
                continue

            # compute node index if necessary
            if not e[0] in nodemap:
                nodemap[e[0]] = idx
                idx += 1

            u = nodemap[e[0]]
            edge = addEdgeLP(pricer.data, u, item)
            v = edge[1]
            assert edge[0] == u
            assert not e[1] in nodemap
            nodemap[e[1]] = v
            node2type[v] = item

            # compute potential values
            assert u in maxp
            maxp[v] = max(maxp[v], maxp[u] + xval * h[item]) if v in maxp else maxp[u] + xval * h[item]

            # count how often item t is packed on node u
            if (u,item) in itemOnNode:
                itemOnNode[u,item] += 1
            else:
                itemOnNode[u,item] = 1

            addEdgeSCIP(model, pricer.data, u, v, False)

        # create the frontier
        for u in [nodemap[i] for i in list(nodemap)]:
            assert u in node2type
            for t in N:
                if (t,node2type[u]) in alpha:
                    number = alpha[t,node2type[u]] if not (u,t) in itemOnNode else alpha[t,node2type[u]] - itemOnNode[u,t]
                    assert number >= 0
                    for k in range(number):
                        e = addEdgeLP(pricer.data, u, t)
                        pricer.data['F'].append(e)

        # compute initial solution
        sol = model.createSol(None)
        for (e,item,xval,yval) in soledgemap:
            if float(yval) < 0.5:
                continue
            u = nodemap[e[0]]
            v = nodemap[e[1]]
            model.setSolVal(sol, pricer.data['var_x'][u,v], xval)
            model.setSolVal(sol, pricer.data['var_y'][u,v], yval)
            model.setSolVal(sol, pricer.data['var_p'][v], maxp[v])

        # add the solution to SCIP
        accepted = model.addSol(sol)
        assert accepted
        model.freeSol(sol)


def runPricing(data, timelim, memlim, lpiter, soledgemap, display, quite):
    """
    solves the semifluid instance with a column generation algorithm on the 'compact' formulation
    output:
    edgemap     - list containing all edge tuples (e,type,idx,x,y)
    stats       - dict containing information about 'status', 'time', 'nnodes', 'nedges', 'primal', 'dual', 'dualroot', 'nxused', 'nyused'
    """
    model  = Model('semi-fluid')
    pricer = fluidPricer()

    createModel(model, pricer, data, soledgemap, lpiter)
    model.includePricer(pricer, 'semi-fluid packing pricer', 'Pricer to compute new important arcs in semi-fluid packing')

    stats = {}
    edgemap = []

    # initialize status
    stats['status']   = 'unknown'
    stats['time']     = 0.0
    stats['nnodes']   = 0
    stats['nedges']   = 0
    stats['nxused']   = 0
    stats['nyused']   = 0
    stats['primal']   = model.infinity()
    stats['dual']     = -model.infinity()
    stats['dualroot'] = -model.infinity()
    stats['gap']      = model.infinity()

    # hide output
    if quite:
        model.hideOutput()

    # set working limits
    model.setPresolve(SCIP_PARAMSETTING.OFF)
    model.setSeparating(SCIP_PARAMSETTING.OFF)
    model.setHeuristics(SCIP_PARAMSETTING.AGGRESSIVE)
    model.setRealParam('limits/time', timelim)
    model.setRealParam('limits/memory', memlim)

    model.optimize()

    if not quite:
        model.printBestSol()

    stats['status']   = model.getStatus()
    stats['time']     = model.getSolvingTime()
    stats['nnodes']   = model.getNNodes()
    stats['dual']     = model.getDualbound()
    stats['dualroot'] = model.getDualboundRoot()
    stats['primal']   = model.getPrimalbound()
    stats['gap']      = model.getGap()

    stats['nedges']   = len(list(pricer.data['var_x']))
    stats['nxused']   = 0
    stats['nyused']   = 0

    # compute edgemap
    for e in list(pricer.data['var_x']):
        xval = model.getVal(pricer.data['var_x'][e])
        yval = model.getVal(pricer.data['var_y'][e])
        t    = pricer.data['LP_edge2item'][e]

        edgemap.append((e, t, xval, yval))

        if xval > 1e-06:
            stats['nxused'] += 1
        if yval > 1e-06:
            stats['nyused'] += 1
    edgemap.sort(key = lambda tup : tup[0][0])

    # display solution
    if display is True:
        displaySolution(data, edgemap)

    # free model and LP
    pricer.data['LP'].free()
    model.free()

    return edgemap, stats


# MAIN
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('usage: python semifluid.py <instance>')
        exit(0)

    data = read_data(str(sys.argv[1]))
    edgemap, stats = runPricing(data, 3600, 10000, -1, None, True, False)

    # print statistics
    print( "%10s | %10s %10s %10s %5s | %10s %10s | %10s %10s %10s" %
           ("status", "primal", "dual", "dualroot", "gap", "time", "nnodes", "nedges", "nxused", "nyused"))
    print( "%10s | %10f %10f %10f %5.1f | %10s %10s | %10s %10s %10s" %
           (stats['status'], \
           stats['primal'], stats['dual'], stats['dualroot'], 100 * stats['gap'], \
           stats['time'], stats['nnodes'], \
           stats['nedges'], stats['nxused'], stats['nyused']))

    edgemap, stats = runPricing(data, 3600, 10000, -1, edgemap, True, False)
    print( "%10s | %10f %10f %10f %5.1f | %10s %10s | %10s %10s %10s" %
           (stats['status'], \
           stats['primal'], stats['dual'], stats['dualroot'], 100 * stats['gap'], \
           stats['time'], stats['nnodes'], \
           stats['nedges'], stats['nxused'], stats['nyused']))
