from pyscipopt import Model, Heur, LP, Conshdlr, quicksum, SCIP_PARAMSETTING, SCIP_RESULT, SCIP_PRESOLTIMING, SCIP_HEURTIMING, SCIP_PROPTIMING
from pyscipopt.linexpr import LinExpr
from data import read_data
from math import floor
import timeit
import sys
import math
import random
class MyHeur(Heur):

    def heurexec(self, heurtiming, nodeinfeasible):
        E_list = self.model.data['E_list']
        E = self.model.data['E']
        nodes = self.model.data['nodes']
        l = self.model.data['l']
        w = self.model.data['w']
        h = self.model.data['h']
        N = self.model.data['N']
        nodetype = self.model.data['nodetype']
        y = self.model.data['y']
        phi = self.model.data['phi']

        score = {}
        for (i,j) in E_list:
            score[i,j] = (10 + l[nodetype[j]])**2 * (w[nodetype[j]] / h[nodetype[j]]) * (self.model.getSolVal(None, y[i,j]))

        E_sorted = sorted(E_list, key=lambda e: score[e], reverse=True)
        taken = {0 : True}
        left = { 0 : self.model.data['D']}
        nused = [0] * len(N)

        for i in nodes:
            left[i] = l[nodetype[i]]

        sol = self.model.createSol(self)

        starts = [0]
        added = True

        # add all fixed variables
        for (i,j) in [(i,j) for (i,j) in E if y[i,j].getLbLocal() > 0.5]:
            left[i] -= l[nodetype[j]]
            taken[j] = True
            nused[nodetype[j]] += 1
            self.model.setSolVal(sol, y[i,j], 1)
            if left[i] < 0:
                return {"result": SCIP_RESULT.DIDNOTFIND}

        # add candidates which do are not used so far
        nrounds = 0
        while len(starts) > 0 and nrounds < 50:
            nrounds += 1
            while len(starts) > 0:
                s = starts.pop()
                candidates = [(i,j) for (i,j) in E_sorted if i == s and abs(y[i,j].getLbLocal() - y[i,j].getUbLocal()) > self.model.feastol()]

                for (i,j) in candidates:
                    if j in taken or l[nodetype[j]] > left[i] or nused[nodetype[j]] > 0:
                        continue
                    taken[j] = True
                    left[i] -= l[nodetype[j]]
                    nused[nodetype[j]] += 1
                    assert left[i] >= 0
                    starts.append(j)
                    self.model.setSolVal(sol, y[i,j], 1)

            # fill the container up
            for (i,j) in [(i,j) for (i,j) in E_sorted if i in taken and not j in taken and abs(y[i,j].getLbLocal() - y[i,j].getUbLocal()) > self.model.feastol()]:
                if left[i] >= l[nodetype[j]] and not j in taken:
                    taken[j] = True
                    left[i] -= l[nodetype[j]]
                    nused[nodetype[j]] += 1
                    self.model.setSolVal(sol, y[i,j], 1)
                    starts.append(j)

        # set solution value of phi to something feasible; this will be adjusted during the CONSCHECK call
        self.model.setSolVal(sol, phi, 0.5 * phi.getLbLocal() + 0.5 * phi.getUbLocal())
        accepted = self.model.trySol(sol, printreason = False)


        if accepted:
            return {"result": SCIP_RESULT.FOUNDSOL}
        else:
            return {"result": SCIP_RESULT.DIDNOTFIND}


class MyConshdlr(Conshdlr):

    def create_lp(self, modeldata):
        # get instance specific data
        s, nnodes, nodetype= modeldata['s'], modeldata['nnodes'], modeldata['nodetype']
        E, H, N, w, h = modeldata['E'], modeldata['H'], modeldata['N'], modeldata['w'], modeldata['h']
        nodes = range(1,nnodes)

        lp = LP()
        x = {}
        p = {}

        # create variables
        lp.addCols([[]] * nnodes, None, [0] * nnodes) # x variables
        lp.addCols([[]] * nnodes, None, [0] * nnodes) # p variables
        assert(lp.ncols() == 2*nnodes)
        for i in nodes:
            x[i] = i - 1
            lp.chgObj(x[i], -w[nodetype[i]])
            p[i] = nnodes + i - 1

        # create and store linear rows
        self.data['cons_enable'] = {}
        self.data['cons_demand'] = {}
        self.data['cons_potential'] = {}
        self.data['cons_potentialroot'] = {}
        self.data['cons_ubx'] = {}
        self.data['cons_ubp'] = {}

        for j in nodes:
            lp.addRow([(x[j], -1)], lhs = 0)
            self.data['cons_enable'][j] = lp.nrows() - 1

        for n in N:
            lp.addRow([(x[i],-1) for i in nodes if nodetype[i] == n], lhs = -1)
            self.data['cons_demand'][n] = lp.nrows() - 1

        for (i,j) in [(i,j) for (i,j) in E if i != s]:
            lp.addRow([(p[j],1), (p[i],-1), (x[j], -h[nodetype[j]])], lhs = -H)
            self.data['cons_potential'][i,j] = lp.nrows() - 1

        for (i,j) in [(i,j) for (i,j) in E if i == s]:
            lp.addRow([(p[j],1), (x[j], -h[nodetype[j]])], lhs = -H)
            self.data['cons_potentialroot'][j] = lp.nrows() - 1

        for i in nodes:
            lp.addRow([(x[i], -1)], lhs = -1)
            self.data['cons_ubx'][i] = lp.nrows() - 1

        for i in nodes:
            lp.addRow([(p[i], -1)], lhs = -H)
            self.data['cons_ubp'][i] = lp.nrows() - 1

        # store constraint handler specific data
        self.data['lp'] = lp
        self.data['x'] = x
        self.data['p'] = p

        # for debugging
        self.nlps = 0

    def update_lp(self, yvals):
        modeldata = self.model.data
        E = modeldata['E']
        lp = self.data['lp']
        nodes = range(1,modeldata['nnodes'])

        # update constraint sides
        for j in nodes:
            lhs = sum(-yvals[i,k] for (i,k) in E if k == j)
            lp.chgSide(self.data['cons_enable'][j], lhs, lp.infinity())
            # print('chg lhs of cons_enable[',j,'] to ', lhs)

        for (i,j) in [(i,j) for (i,j) in E if i != modeldata['s']]:
            lp.chgSide(self.data['cons_potential'][i,j], (yvals[i,j] - 1) * modeldata['H'], lp.infinity())
            # print('chg lhs of cons_potential[',i,j,'] to ', (yvals[i,j] - 1) * modeldata['H'])

        for (i,j) in [(i,j) for (i,j) in E if i == modeldata['s']]:
            lp.chgSide(self.data['cons_potentialroot'][j], (yvals[i,j] - 1) * modeldata['H'], lp.infinity())
            # print('chg lhs of cons_potentialroot[',j,'] to ', (yvals[i,j] - 1) * modeldata['H'])

    def consinit(self, constraints):
        # create LP
        self.data = {}
        self.create_lp(self.model.data)

        self.counter = 0
        self.y0 = {}
        for (i,j) in self.model.data['E']:
            self.y0[i,j] = 1

        assert('lp' in self.data)
        assert('x' in self.data)
        assert('p' in self.data)


    def consexit(self, constraints):
        # free LP
        del self.data['lp']
        self.data = {}

    def conscheck(self, constraints, solution, checkintegrality, checklprows, printreason):
        # print("[conscheck]")
        modeldata = self.model.data
        lp = self.data['lp']

        # collect y values
        y = modeldata['y']
        yvals = {}
        for (i,j) in modeldata['E']:
            yvals[i,j] = self.model.getSolVal(solution, modeldata['y'][i,j])

        # update LP
        self.update_lp(yvals)
        objval = lp.solve()

        assert lp.isPrimalFeasible()
        assert lp.isDualFeasible()

        primalsol = lp.getPrimal()
        xvals = {}
        pvals = {}

        #     # violation of the solution arises from a too small z* value; adjust z to the solution value of the LP to obtain
        #     # a primal feasible solution
        #     # @todo check whether heuristics are fine with the change of the solution value
        if objval < self.model.getPrimalbound() - self.model.feastol():
            self.model.setSolVal(solution, modeldata['phi'], objval)
            return {"result": SCIP_RESULT.FEASIBLE}
        return {"result": SCIP_RESULT.INFEASIBLE}

    def compute_cut(self, yvals, dualsol, objval, checkassert):
        modeldata = self.model.data
        E = modeldata['E']
        H = modeldata['H']
        y = modeldata['y']

        cut = quicksum([])
        activity = 0

        for j in modeldata['nodes']:
            assert dualsol[self.data['cons_enable'][j]] >= -self.model.feastol()
            cut += quicksum(-y[i,k] * dualsol[self.data['cons_enable'][j]] for (i,k) in E if k == j)
            activity += sum(-yvals[i,k] * dualsol[self.data['cons_enable'][j]] for (i,k) in E if k == j)

        # compute a Benders cut
        for n in modeldata['N']:
            assert dualsol[self.data['cons_demand'][n]] >= -self.model.feastol()
        cut += quicksum(-1 * dualsol[self.data['cons_demand'][n]] for n in modeldata['N'])
        activity += sum(-1 * dualsol[self.data['cons_demand'][n]] for n in modeldata['N'])

        for (i,j) in [(i,j) for (i,j) in E if i != modeldata['s']]:
            assert dualsol[self.data['cons_potential'][i,j]] >= -self.model.feastol()
        cut += quicksum(dualsol[self.data['cons_potential'][i,j]] * (y[i,j] - 1) * H for (i,j) in E if i != modeldata['s'])
        activity += sum(dualsol[self.data['cons_potential'][i,j]] * (yvals[i,j] - 1) * H for (i,j) in E if i != modeldata['s'])

        for (i,j) in [(i,j) for (i,j) in E if i == modeldata['s']]:
            assert dualsol[self.data['cons_potentialroot'][j]] >= -self.model.feastol()
        cut += quicksum(dualsol[self.data['cons_potentialroot'][j]] * (y[i,j] - 1) * H for (i,j) in E if i == modeldata['s'])
        activity += sum(dualsol[self.data['cons_potentialroot'][j]] * (yvals[i,j] - 1) * H for (i,j) in E if i == modeldata['s'])

        cut += quicksum(-1 * dualsol[self.data['cons_ubx'][i]] for i in modeldata['nodes'])
        activity += sum(-1 * dualsol[self.data['cons_ubx'][i]] for i in modeldata['nodes'])

        cut += quicksum(-H * dualsol[ self.data['cons_ubp'][i]] for i in modeldata['nodes'])
        activity += sum(-H * dualsol[self.data['cons_ubp'][i]] for i in modeldata['nodes'])

        if checkassert:
            assert(abs(activity - objval) <= self.model.feastol())

        if activity > self.model.getSolVal(None, modeldata['phi']) + self.model.feastol():
            self.model.addCons(cut - modeldata['phi'] <= 0, modifiable = True, local = True, name='cut_' + str(self.counter))
            # self.model.writeTransformedProblem()
            self.counter += 1

            return SCIP_RESULT.CONSADDED
        else:
            return SCIP_RESULT.FEASIBLE


    def consenfolp(self, constraints, nusefulconss, solinfeasible):
        # print("[consenfolp]")
        modeldata = self.model.data
        lp = self.data['lp']

        # collect y values
        y = modeldata['y']
        yvals = {}
        for (i,j) in modeldata['E']:
            yvals[i,j] = self.model.getSolVal(None, modeldata['y'][i,j])

        # update LP
        self.update_lp(yvals)
        objval = lp.solve()
        assert lp.isPrimalFeasible()
        assert lp.isDualFeasible()

        # check whether the LP solution is feasible
        if self.model.getSolVal(None, modeldata['phi']) >= objval - self.model.feastol():
            return {"result": SCIP_RESULT.FEASIBLE}

        nodes = modeldata['nodes']
        dualsol = lp.getDual()


        # compute a Benders cut
        result = self.compute_cut(yvals, dualsol, objval, True)

        # stop if the current solution is feasible
        if result == SCIP_RESULT.FEASIBLE:
            return {"result": SCIP_RESULT.FEASIBLE}

        # # compute cut between y* and y0
        # for (i,j) in modeldata['E']:
        #     self.y0[i,j] = 0.5 * yvals[i,j] + 0.5 * self.y0[i,j]
        # self.update_lp(self.y0)
        # result = self.compute_cut(self.y0, dualsol, objval, False)
        # # print(result == SCIP_RESULT.FEASIBLE)

        # # compute cut between y* and the best solution so far
        # if self.model.getNSols() > 0:
        #     bestsol = self.model.getBestSol()
        #     for (i,j) in modeldata['E']:
        #         self.y0[i,j] = 0.5 * yvals[i,j] + 0.5 * self.model.getSolVal(bestsol, modeldata['y'][i,j])
        #     self.update_lp(self.y0)
        #     result = self.compute_cut(self.y0, dualsol, objval, False)
        #     # print(result == SCIP_RESULT.FEASIBLE)

        return {"result": SCIP_RESULT.CONSADDED}


    def consenfops(self, constraints, nusefulconss, solinfeasible, objinfeasible):
        if self.model.getSolvingTime() > self.model.data['timelim']:
            return {"result": SCIP_RESULT.DIDNOTFIND}
        return self.consenfolp(constraints, nusefulconss, solinfeasible)


    def conslock(self, constraint, nlockspos, nlocksneg):
        # lock all y variables in both directions
        y = self.model.data['y']
        phi = self.model.data['phi']
        for (i,j) in list(y):
            self.model.addVarLocks(y[i,j], nlocksneg + nlockspos, nlocksneg + nlockspos)

        # phi just needs a down lock
        self.model.addVarLocks(phi, nlocksneg + nlockspos, 0)


def create_master(data, s, nnodes, nodetype, E):
    model = Model('benders')
    nodes = range(1,nnodes)
    l = data['l']
    y = {}

    # create conshdlr and include it to SCIP
    conshdlr = MyConshdlr()
    model.includeConshdlr(conshdlr, "PyCons", "custom constraint handler implemented in python",
                      sepapriority = 1, enfopriority = 100, chckpriority = -1000000, sepafreq = 1, propfreq = -1,
                      eagerfreq = 1, maxprerounds = -1, delaysepa = False, delayprop = False, needscons = False,
                      presoltiming = SCIP_PRESOLTIMING.FAST, proptiming = SCIP_PROPTIMING.BEFORELP)

    # create variables
    for (i,j) in E:
        y[i,j] = model.addVar('y_' + str(i) + ',' + str(j), obj=0, lb=0, ub=1, vtype='B')

    phi = model.addVar('phi', obj=1.0, lb=-model.infinity(), ub=model.infinity(), vtype='C')

    model.data = data
    model.data['s'] = s
    model.data['nnodes'] = nnodes
    model.data['nodes'] = range(1,nnodes)
    model.data['nodetype'] = nodetype
    model.data['E'] = E
    model.data['E_list'] = sorted(list(E))
    model.data['y'] = y
    model.data['phi'] = phi

    # sum_{e in delta^{-}(i)} y_e <= 1
    for j in nodes:
        model.addCons(quicksum(y[i,j] for (i,k) in E if k == j) <= 1, name='incoming_' + str(j))

    # sum_{e in delta^{+}(s) l[e] * y_e <= D}
    model.addCons(quicksum(l[nodetype[i]] * y[s,i] for i in nodes) <= data['D'], name='depth')

    # sum_{e in delta^{+}(i)} l(e) * y_e <= sum_{e in delta^{-}(i)} l(e) * y_e
    for i in nodes:

        # constraint is redundant for nodes without outgoing edges
        if nodetype[i] == data['N'][0]:
            continue

        sum1 = quicksum(l[nodetype[j]] * y[i,j] for (k,j) in E if k == i)
        sum2 = quicksum(-l[nodetype[i]] * y[j,i] for (j,k) in E if k == i)
        model.addCons(sum1 + sum2 <= 0, name='length_' + str(i))

    # include the heuristic
    heuristic = MyHeur()
    model.includeHeur(heuristic, "PyHeur", "custom heuristic implemented in python", "B", freq=1, timingmask=SCIP_HEURTIMING.DURINGLPLOOP)

    return model


def runBenders(data, timelim, memlim, display, quite):
    """
    creates and solves instance via Bender's Decomposition
    output:
    edgemap     - list containing all edge tuples (e,type,x,y)
    stats       - dict containing information about 'status', 'time', 'nnodes', 'readtime', 'nedges', 'primal', 'dual', 'dualroot', 'nxused', 'nyused'
    """
    D, W, H, N, l, h, w = data['D'], data['W'], data['H'], data['N'], data['l'], data['h'], data['w']
    stats = {}

    # create graph
    s, nnodes, nodetype, E = create_graph(N, D, l)

    # create master problem
    model = create_master(data, s, nnodes, nodetype, E)

    # set working limits
    model.setRealParam('limits/time', timelim)
    model.setRealParam('limits/memory', memlim)
    model.setIntParam('display/freq', 10000)
    # model.setHeuristics(SCIP_PARAMSETTING.AGGRESSIVE)

    # store time limit for proper termination
    model.data['timelim'] = timelim

    # hide output
    if quite:
        model.hideOutput()

    model.optimize()
    if not quite:
        model.printBestSol()
        model.printStatistics()

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

    # if display:

    return stats


def create_graph(N, D, l):
    V     = {}
    nodetype = {}
    E     = []
    s = 0
    nnodes = 1

    # create all nodes
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

    return s, nnodes, nodetype, E

# MAIN
if __name__ == '__main__':
    if len(sys.argv) < 4:
        print('usage: python benders.py <instance> <time limit> <memory limit> [display] [quite]')
        exit(0)

    # read input
    data = read_data(str(sys.argv[1]))
    timelim = float(sys.argv[2])
    memlim = float(sys.argv[3])

    # optional parameter
    display = True if len(sys.argv) >= 5 and str(sys.argv[4]) == 'True' else False
    quite = True if len(sys.argv) >= 6 and str(sys.argv[5]) == 'True' else False

    # solve problem
    stats = runBenders(data, timelim, memlim, display, quite)
    print('STATS: ' + str(stats))
