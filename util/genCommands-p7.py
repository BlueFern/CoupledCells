import argparse
import math

# The output from this script can be redirected to a file.

if __name__ == "__main__":

    argParser = argparse.ArgumentParser(
        description='Generate output containing N commands to be used by batecher with -np N where N is the number of procs.')
    argParser.add_argument('numProcs', type=int, help='Number of processors to be used with batcher.')
    argParser.add_argument('numSteps', type=int, help='Number of steps to be divided between the specified number of processors.')
    argParser.add_argument('circQuads', type=int)
    argParser.add_argument('axialQuads', type=int)
    argParser.add_argument('circScale', type=int)
    argParser.add_argument('axialScale', type=int)
    args = argParser.parse_args()

    # The ceiling value for the number of steps per processor.
    stepPerProc = int(math.ceil(args.numSteps / float(args.numProcs)))

    # For every processor calculate the start and end steps.
    for proc in range(args.numProcs):
        start = stepPerProc * proc
        adjustedStep = stepPerProc

        # Make sure the last step is set to not overstep the total number of steps.
        if proc == args.numProcs - 1:
            adjustedStep = args.numSteps - start

        end = start + adjustedStep - 1

        # Adjust relative path as necessary.
        print 'python2.7 ../../CoupledCells/util/hdf5ToVTU.py ' + str(start) + ' ' + str(end) + ' ' + str(args.circQuads) + ' ' + str(args.axialQuads) + ' ' + str(args.circScale) + ' ' + str(args.axialScale) + " smc"
        print 'python2.7 ../../CoupledCells/util/hdf5ToVTU.py ' + str(start) + ' ' + str(end) + ' ' + str(args.circQuads) + ' ' + str(args.axialQuads) + ' ' + str(args.circScale) + ' ' + str(args.axialScale) + " ec"