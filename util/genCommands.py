import argparse
import math

if __name__ == "__main__":

    argParser = argparse.ArgumentParser(
        description='Generate output containing N commands to be used by batecher with -np N where N is the number of procs.')
    argParser.add_argument('numProcs', type=int, help='Number of processors to be used with batcher.')
    argParser.add_argument('numSteps', type=int, help='Number of steps to be divided between the specified number of processors.')
    args = argParser.parse_args()

    # The ceiling value for the number of steps per processor.
    stepPerProc = int(math.ceil(args.numSteps / float(args.numProcs)))

    for proc in range(args.numProcs):
        start = stepPerProc * proc
        adjustedStep = stepPerProc
        if proc == args.numProcs - 1:
            adjustedStep = args.numSteps - start
        end = start + adjustedStep - 1

        print 'time python2.7 ../scripts/bifurcation_ec_hdf5ToVTU.py ' + str(start) + ' ' + str(end)
        print 'time python2.7 ../scripts/bifurcation_smc_hdf5ToVTU.py ' + str(start) + ' ' + str(end)