from data import read_data
from semifluid import runPricing
# from compact import runCompact
# from compact2 import runCompact2
from compact_improve import runCompactImprove
from benders import runBenders
from evalinstance import computeNEdgesComplex
from math import floor
from os import listdir
from os.path import isfile, join
import sys


def GE(x,y):
    return float(x) - float(y) >= -1e-06


def printStatistics(stats):
    print( "%10s | %10f %10f %10f %5.1f | %10s %10s | %10s %10s %10s" %
           (stats['status'], \
            stats['primal'], stats['dual'], stats['dualroot'], 100 * stats['gap'], \
            stats['time'], stats['nnodes'], \
            stats['nedges'], stats['nxused'], stats['nyused']))


def checkStatistic(stats1, stats2):
    if not GE(stats1['primal'], stats2['dual']) or not GE(stats2['primal'], stats1['dual']):
        print("ERROR: primalbound < dualbound")


def runInstance(instance):
    readtimelim = 120
    timelim     = 60.0
    memlim      = 10000
    lpiter      = -1
    quite       = True

    data = read_data(instance)

    # # run compact1 model
    # compact_edgemap, compact_stats = runCompact(data, timelim, memlim, readtimelim, None, "compact.sol", False, quite)
    # printStatistics(compact_stats)

    # # run compact2 model
    # compact2_stats  = runCompact2(data, timelim, memlim, False, quite)
    # printStatistics(compact2_stats)

    # checkStatistic(compact_stats, compact2_stats)

    stats_compact = runCompactImprove(data, timelim, memlim, False, quite)
    printStatistics(stats_compact)
    stats_benders = runBenders(data, timelim, memlim, False, quite)
    printStatistics(stats_benders)
    checkStatistic(stats_compact, stats_benders)


    # # run compact model (+ warmstart)
    # compact_edgemap, compact_stats = runCompact(data, timelim, memlim, readtimelim, None, "compact.sol", False, quite)
    # printStatistics(compact_stats)
    # compact_edgemap, compact_stats = runCompact(data, timelim, memlim, readtimelim, "compact.sol", None, False, quite)
    # printStatistics(compact_stats)

    # # run pricing model (+warmstart)
    # pricer_edgemap, pricer_stats  = runPricing(data, timelim, memlim, lpiter, None, False, quite)
    # printStatistics(pricer_stats)
    # checkStatistic(compact_stats, pricer_stats)

    # if float(pricer_stats['primal']) <= float(compact_stats['primal']):
    #     pricer_edgemap, pricer_stats  = runPricing(data, timelim, memlim, lpiter, pricer_edgemap, False, quite)
    # else:
    #     pricer_edgemap, pricer_stats  = runPricing(data, timelim, memlim, lpiter, compact_edgemap, False, quite)
    # printStatistics(pricer_stats)


if __name__ == '__main__':
    if len(sys.argv) <= 1:
        print("usage: python run.py {easy,hard,<instance>}")
        exit(0)

    # print head line
    print( "%10s | %10s %10s %10s %5s | %10s %10s | %10s %10s %10s" %
           ("status", "primal", "dual", "dualroot", "gap", "time", "nnodes", "nedges", "nxused", "nyused"))

    if str(sys.argv[1]) == "easy" or str(sys.argv[1]) == "hard":
        reader = open(str(sys.argv[1]) + ".test", 'r')
        for line in reader:
            instance = str(line.split()[0])
            nedges   = float(line.split()[1])
            print("%s %e" % (instance, nedges))
            runInstance("data/" + instance)
        reader.close()
    else:
        instance = str(sys.argv[1])
        data = read_data(instance)
        runInstance(instance)

