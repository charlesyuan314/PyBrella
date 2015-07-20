#!/bin/env python
# pybrella.py
# PyBrella by Charles Yuan, 2014
# PyBrella is a script designed to automate the umbrella sampling process, using the AMBER molecular dynamics package.

# PyBrella automates the setup and execution of MD runs for umbrella sampling on the distance parameter, using user-
# defined start and end points. After completion of the initial sampling, PyBrella performs confirmation of sample
# convergence based on the 80/20 chi-squared distance method. It also examines the sample histograms and automatically
# begins new samples where necessary to reach a specified sampling threshold.

# Contains modified code from the pdbTools project by Mike Harms: https://code.google.com/p/pdb-tools/

from __future__ import division, print_function
from argparse import ArgumentParser
from atom_tools import atomDist, calcCenterAtoms
from glob import glob
from histogram import histogram
from shared import *
from stat import S_IEXEC
from time import sleep, strftime

__author__ = "Charles Yuan"
__credits__ = ["Charles Yuan", "David Koes", "Matthew Baumgartner"]
__license__ = "GPL"
__version__ = "3.0"
__maintainer__ = "Charles Yuan"
__email__ = "charlesyuan314@gmail.com"
__status__ = "Development"


def log(message, toFile=True):
    """Log a message, optionally to the args.log file."""
    sys.stderr.write(message)
    if toFile:
        try:
            if args.log is not None:
                ANSIlist = [BOLD, UNDERLINE, BLACK, RED, GREEN, YELLOW, BLUE, MAGENTA, CYAN, WHITE, END]
                newMessage = message
                for code in ANSIlist:
                    newMessage = newMessage.replace(code, '')
                with open(args.log, 'a') as logFile:
                    logFile.write(newMessage)
        except NameError:
            pass


def main():
    """Calls major subroutines, prepares files, and prints messages. Hands control over to eventLoop afterwards."""
    parse()

    runInit = not prep()

    log(CYAN + BOLD + "PyBrella" + END + CYAN + " version " + MAGENTA + "%s" % __version__ + CYAN + " starting up on "
        + MAGENTA + "%s\n" % args.queue_name + END)
    log("Configured for " + MAGENTA + "%i" % args.precision + END + " decimal places.\n" + END)
    log("Using topology at " + UNDERLINE + "%s\n" % PRMTOPPATH + END)
    log("Using coordinate folder at " + UNDERLINE + "%s\n" % COORDPATH + END)
    log("Using distance file at " + UNDERLINE + "%s\n" % DISTPATH + END)
    log("Using " + MAGENTA + "%i" % args.steps + END + " timesteps per simulation.\n")
    log(BLUE + "Simulation will occur over " + MAGENTA + "%.3f" % MAXDIST + BLUE +
        " angstroms, incrementing every " + MAGENTA + "%.3f" % INCRDIST + BLUE + " angstroms.\n" + END)
    log("Burst mode is set to " + MAGENTA + "%i" % args.burst + END + " runs per set.\n")
    log("Chi-squared analysis is ")
    if args.chi_squared > 0:
        log(MAGENTA + "enabled" + END + " at " + MAGENTA + "%.*f" % (args.precision, args.chi_squared) + END + ".\n")
    else:
        log(MAGENTA + "disabled" + END + ".\n")
    log("Normalization is ")
    if args.normalize:
        log(MAGENTA + "enabled" + END + ".\n")
    else:
        log(MAGENTA + "disabled" + END + ".\n")
    log("Other parameters: force = %.*f; bin = %.3f; thres = %i; harmonic = %.3f\n" %
        (args.precision, args.force, BINWIDTH, args.thres, HARMONIC))

    determineRestraint()  # Determine which atoms will be used for restaints

    if runInit:
        init()
    else:
        finishedRuns = checkFinishedRuns()
        if finishedRuns:
            # Check chi-squared
            if args.chi_squared > 0:
                chiSquaredPassed = checkChiSquared()
                if chiSquaredPassed:
                    analysis()
            else:
                analysis()
        else:
            log(CYAN + "Runs have not yet finished.\n" + END)
    if args.loop > 0:
        eventLoop()


def parse():
    """Prepare the argument parser."""
    parser = ArgumentParser(description="execute umbrella sampling with AMBER.")
    parser.add_argument("--prmtop", '-p', help="AMBER topology file", type=str, action=FullPath, required=True)
    parser.add_argument("--coord", '-c', help="coordinate folder of entire trajectory", type=str, action=FullPath,
                        default="trajectory")
    parser.add_argument("--dist", '-d', help="dist_vs_t file, sorted by distance", type=str, action=FullPath,
                        default="sorted_dist_vs_t")
    parser.add_argument("--steps", help="nstlim value for AMBER config (default is 50000)", type=int,
                        action="store", default=50000)

    parser.add_argument("--init", help="force the initialization process for adding new runs.", action="store_true")

    parser.add_argument("--max", help="magnitude, in angstroms, of distance the ligand must travel. Must be a multiple"
                                      " of the increment.",
                        action="store", type=float, required=True)
    parser.add_argument("--incr", help="positive length increment for initial runs in angstroms (default is 1)",
                        action="store", type=float, default=1)

    parser.add_argument("--force", '-f', help="force constant (default is 5)", action="store", type=float, default=5)
    parser.add_argument("--bin",
                        help="bin width for histograms in analysis runs. >=0.01 recommended.",
                        action="store", type=float, default=0.01)
    parser.add_argument("--thres", '-t', help="threshold for histogram sampling during analysis", action="store",
                        type=int, default=5000)
    parser.add_argument("--chi_squared",
                        help="for analysis, continue a run if the last 20%% deviates from the first 80%% by this "
                             "distance. Set to 0 to disable checking. (default is 0.05)",
                        type=float, action="store", default=0.05)
    parser.add_argument("--harmonic", help="left and right widths of harmonic width (default is 15)", type=float,
                        default=15)
    parser.add_argument("--burst", '-b', help="number of runs per burst (default is 10)", type=int, default=10)
    parser.add_argument("--normalize", '-n', help="ensure that all run sets have the same number of burst runs.",
                        action="store_true")

    parser.add_argument("--queue_name", '-q', help="queue name", type=str, required=True)
    parser.add_argument("--precision", help="number of decimal places used in calculations. "
                                            "Note that distances are limited by AMBER to 3 places.",
                        type=int, action="store", default=6)
    parser.add_argument("--log", help="log output file", type=str, action=FullPath)
    parser.add_argument("--hist", help="specify this path to output the analysis histogram.", type=str,
                        action=FullPath)
    parser.add_argument("--loop", '-l', help="number of seconds before the loop checks status. "
                                             "Set to 0 to disable the loop.", type=int, action="store", default=60)
    global args
    args = parser.parse_args()


def prep():
    """Sets up variables for the entire script. Also detects whether analysis runs are necessary and
    returns a boolean to indicate this."""
    global WORKDIR
    WORKDIR = os.getcwd()

    global PAUSE
    if args.loop <= 0:
        PAUSE = 1
    else:
        PAUSE = args.loop

    global STARTDIST
    STARTDIST = 0

    global AVERAGEHILLWIDTH
    AVERAGEHILLWIDTH = 0

    global PRMTOPPATH
    global COORDPATH
    global DISTPATH
    global HISTPATH

    # Path modifications
    if args.prmtop is not None and not os.path.isabs(args.prmtop):
        PRMTOPPATH = WORKDIR + "/" + args.prmtop
    else:
        PRMTOPPATH = args.prmtop
    if not os.path.isfile(PRMTOPPATH):
        fail(RED + UNDERLINE + "Error:" + END + RED + " specified topology file does not exist.\n" + END)
    if args.coord is not None and not os.path.isabs(args.coord):
        COORDPATH = WORKDIR + "/" + args.coord
    else:
        COORDPATH = args.coord
    if not os.path.isdir(COORDPATH):
        fail(RED + UNDERLINE + "Error:" + END + RED + " specified coordinate path is not a directory.\n" + END)
    if args.hist is not None and not os.path.isabs(args.hist):
        HISTPATH = WORKDIR + "/" + args.hist
    else:
        HISTPATH = args.hist
    if args.dist is not None and not os.path.isabs(args.dist):
        DISTPATH = WORKDIR + "/" + args.dist
    else:
        DISTPATH = args.dist
    if not os.path.isfile(DISTPATH):
        fail(RED + UNDERLINE + "Error:" + END + RED + " specified distance file does not exist.\n" + END)

    # Round AMBER distance parameters to 3 decimal places
    global MAXDIST
    MAXDIST = round(args.max, 3)
    global INCRDIST
    INCRDIST = round(args.incr, 3)
    global BINWIDTH
    BINWIDTH = round(args.bin, 3)
    global HARMONIC
    HARMONIC = round(args.harmonic, 3)

    global RUNANALYSIS
    RUNANALYSIS = False
    if [float(name.split('/')[-1]) for name in glob("%s/*.*" % WORKDIR) if os.path.isdir(os.path.join(WORKDIR, name))]:
        # There's already stuff here
        RUNANALYSIS = True
    if args.init:
        RUNANALYSIS = False

    global PROTEINATOM
    PROTEINATOM = 0
    global LIGANDATOM
    LIGANDATOM = 0

    return RUNANALYSIS


class RunSet(object):
    """A RunSet object holds a burst set of runs. Each run holds
     its own topology, input, output, and coordinate files. It can execute members of itself
     with the execute(IDs) method. It can fetch frames from the coordinate files and prepare files."""

    _dist = 0  # distance of the ligand from the protein
    _restartFile = None  # path to restart file in trajectory folder

    def __init__(self, dist, restartFile=None):
        self._dist = dist
        self._restartFile = restartFile

    def __str__(self):
        return "<RunSet: dist %.3f, restart at %s>" % (self._dist, self._restartFile)

    @property
    def path(self):  # directory holding run folders
        return "%s/%.3f" % (WORKDIR, self._dist)

    def execute(self, IDs=None):
        """Call runs from their scripts"""
        if IDs is None:
            IDs = sorted([name for name in os.listdir(self.path) if os.path.isdir(os.path.join(self.path, name))])
        for ID in IDs:
            with directory(self.path + "/" + str(ID)):
                if not os.path.isfile("./run.sh"):
                    log(RED + UNDERLINE + "Error:" + END + RED + " run.sh does not exist at %s/%s\n" %
                        (self.path, str(ID)) + END)

                system("./run.sh > /dev/null 2> /dev/null")
        if IDs:
            log(GREEN + "Executed runs for distance %.3f on %s.\n" % (self._dist, strftime("%c") + END))
        return

    def create(self, frames, start=0):
        """Creates the folder structure and writes scripts for a run, given the frame info.
        Returns list of new runs and last index reached."""
        global DISTPATH
        finalDist = 0
        lastDist = 0
        newStart = start
        for index, line in enumerate(frames[0:]):
            thisDist = float(line.split()[1])
            if thisDist > self._dist:
                if thisDist - self._dist < self._dist - lastDist:
                    # choose thisDist
                    self._restartFile = "%s/%.3f.rst" % (COORDPATH, thisDist)
                    finalDist = thisDist
                else:
                    # choose lastDist
                    self._restartFile = "%s/%.3f.rst" % (COORDPATH, lastDist)
                    finalDist = lastDist
                newStart = index
                break
            else:
                lastDist = thisDist
                continue
        difference = finalDist - self._dist
        if abs(difference) > 2:
            log(YELLOW + UNDERLINE + "Warning:" + END + YELLOW + " no appropriate frame"
                                                                 " could be found for %.3f.\n" % self._dist + END)
            return [], newStart
        else:
            if abs(difference) > 0.001:
                log("Created runs near %.3f" % self._dist + " using frame at " + MAGENTA + "%.3f" % finalDist + END +
                    " (difference of " + MAGENTA + "%.3f" % difference + END + ").\n")
            else:
                log("Created runs near %.3f" % self._dist + " using frame at " +
                    MAGENTA + "%.3f" % self._dist + END + ".\n")
            self._dist = finalDist  # override the distance
            if os.path.exists(self.path):
                return self.restart(), 0
            else:
                os.mkdir(self.path)
                with directory(self.path):
                    for i in xrange(args.burst):
                        os.mkdir(self.path + "/" + str(i))
                        self.writeScripts(i)
                return range(args.burst), newStart

    def restart(self):
        """Create new folders for restart runs. Returns list of new runs."""
        self._dist = float(self.path.split('/')[-1])
        self._restartFile = "%s/%.3f.rst" % (COORDPATH, self._dist)
        newID = max(sorted([int(name) for name in os.listdir(self.path) if
                            os.path.isdir(os.path.join(self.path, name))])) + 1
        for i in xrange(newID, newID + args.burst):
            os.mkdir(self.path + "/" + str(i))
            self.writeScripts(i)
        return range(newID, newID + args.burst)

    def writeScripts(self, burstNum):
        """Write the qsub and qscript scripts for this run"""
        global PRMTOPPATH
        with directory(self.path + "/" + str(burstNum)):
            if args.queue_name in ["dept_gpu", "any_gpu", "bahar_gpu"]:
                with open("run.sh", 'w') as runScript:
                    runScript.write("#!/bin/bash\n")
                    runScript.write("sleep 0.25; touch executed; qsub -d . -q %s -S /bin/bash -N PB%.3f_%i -l "
                                    "nodes=1:ppn=1:gpus=1 %s/%s/qscript\n" %
                                    (args.queue_name, self._dist, burstNum, self.path, str(burstNum)))
                with open("qscript", 'w') as qscript:
                    qscript.write("""#!/bin/bash
AMBERHOME=/usr/local/amber14
PATH=/usr/local/amber14/bin:$PATH
pmemd.cuda -O -i umb.in -o umb.out -p %s -c %s -r end.rst -x coord.nc -ref %s
touch finished
rm mdinfo
$cmd
""" % (PRMTOPPATH, self._restartFile, self._restartFile))

            elif args.queue_name == "gpu_short":
                with open("run.sh", 'w') as runScript:
                    runScript.write("#!/bin/bash\n")
                    runScript.write("sleep 0.25; touch executed; qsub -d . -q gpu_short -S /bin/bash -N PB%.3f_%i -l "
                                    "nodes=1:ppn=1:gpus=1 -l feature=titan -l walltime=23:59:59 %s/qscript\n" %
                                    (self._dist, burstNum, self.path))

                with open("qscript", 'w') as qscript:
                    qscript.write("""#!/bin/bash
module purge
module load intel/2013.0
module load amber/14-intel-2013-cuda-5.0
pmemd.cuda -O -i umb.in -o umb.out -p %s -c %s -r end.rst -x coord.nc -ref %s
touch finished
rm mdinfo
$cmd
""" % (PRMTOPPATH, self._restartFile, self._restartFile))

            os.chmod("run.sh", os.stat("run.sh").st_mode | S_IEXEC)

            with open("umb.in", 'w') as inputFile:
                inputFile.write("""&cntrl
 imin = 0, ntx = 1, irest = 0,
 ntpr = 10000, ntwr = 100000, ntwx = 5000, ntxo = 1,
 ntf = 2, ntc = 2, cut = 8.0,
 ntb = 2,  nstlim = %i, dt = 0.002,
 temp0 = 300.0, ntt = 3, ig = -1,
 gamma_ln = 1,
 ntp = 1, pres0 = 1.0, taup = 5.0,
 nmropt = 1, ioutfm = 1,
/

&wt TYPE='DUMPFREQ', istep1=1, /
&wt TYPE='END', /
DISANG=dist.rst
DUMPAVE=dist_vs_t
""" % args.steps)
            with open("dist.rst", 'w') as distFile:
                global PROTEINATOM
                global LIGANDATOM
                distFile.write("&rst iat=%i, %i, r1=%.3f, r2=%.3f, r3=%.3f, "
                               "r4=%.3f, rk2=%.*f, rk3=%.*f /\n" %
                               (PROTEINATOM, LIGANDATOM, self._dist - HARMONIC, self._dist, self._dist,
                                self._dist + HARMONIC, args.precision, args.force, args.precision, args.force))


def init():
    """Begin initial runs."""
    global WORKDIR
    global PRMTOPPATH
    global COORDPATH
    global DISTPATH

    log(CYAN + "Beginning initialization.\n" + END)

    # Generate list of runs
    calcInitDist()
    if args.init:  # may already be runs, don't mess with them
        dists = ["%.3f" % dist for dist in frange(STARTDIST, STARTDIST + MAXDIST + INCRDIST, INCRDIST)]
        preexisting = [name.split('/')[-1] for name in glob("%s/*.*" % WORKDIR)
                       if os.path.isdir(os.path.join(WORKDIR, name))]
        dists = [float(dist) for dist in dists if dist not in preexisting]
    else:
        dists = frange(STARTDIST, STARTDIST + MAXDIST + INCRDIST, INCRDIST)
    if dists:
        with open(DISTPATH, 'r') as dist_vs_t:
            frames = dist_vs_t.read().splitlines()
        start = 0
        for dist in dists:
            newRunSet = RunSet(dist)
            runIDs, start = newRunSet.create(frames, start)
            newRunSet.execute(IDs=runIDs)
        log(GREEN + "Initial runs have begun.\n" + END)
    else:
        log(GREEN + "No initial runs necessary.\n" + END)


def determineRestraint():
    """Determine the atom groups used for the center-of-mass restaints and sets them in the PROTEINATOM and LIGANDATOM
    lists."""
    global PROTEINATOM
    global LIGANDATOM
    global PRMTOPPATH
    global COORDPATH

    with open(COORDPATH + "/initial.pdb", 'r') as p:
        LIGANDATOM, PROTEINATOM, ligandAtomCoord, proteinAtomCoord, ligandCenter, proteinCenter = \
            calcCenterAtoms(p.readlines())
    log("Selected ligand center atom ID " + MAGENTA + "%i" % LIGANDATOM + END +
        " at (%.3f, %.3f, %.3f), " % (ligandAtomCoord[0], ligandAtomCoord[1], ligandAtomCoord[2]) +
        "closest to (%.3f, %.3f, %.3f).\n" % (ligandCenter[0], ligandCenter[1], ligandCenter[2]))
    log("Selected protein center atom ID " + MAGENTA + "%i" % PROTEINATOM + END +
        " at (%.3f, %.3f, %.3f), " % (proteinAtomCoord[0], proteinAtomCoord[1], proteinAtomCoord[2]) +
        "closest to (%.3f, %.3f, %.3f).\n" % (proteinCenter[0], proteinCenter[1], proteinCenter[2]))


def calcInitDist():
    """Calculates the "initial distance" based on the first frame of the coordinate file."""
    global RUNANALYSIS
    global PRMTOPPATH
    global COORDPATH
    global STARTDIST
    global PROTEINATOM
    global LIGANDATOM
    if not RUNANALYSIS:
        log("Calculating initial distance.\n")
        with open(COORDPATH + "/initial.pdb", 'r') as pdb:
            STARTDIST = atomDist(pdb.readlines(), LIGANDATOM, PROTEINATOM)
    else:  # analysis
        STARTDIST = min([float(name.split('/')[-1]) for name in glob("%s/*.*" % WORKDIR)
                         if os.path.isdir(os.path.join(WORKDIR, name))])
    if not RUNANALYSIS:
        log(BLUE + "Initial distance is " + MAGENTA + "%.3f" % STARTDIST + BLUE + ". " + END)
        if MAXDIST is not None:
            log(BLUE + "Analysis endpoint will be " + MAGENTA + "%.3f" %
                (MAXDIST + STARTDIST) + BLUE + " angstroms.\n" + END)
        else:
            log('\n')


def checkNormalized():
    """Check if all RunSets have the same number of burst runs"""
    log(BLUE + "Checking if run sets have the same number of burst runs.\n" + END)

    # Get list of subdirectories
    distances = [float(name.split('/')[-1]) for name in glob("%s/*.*" % WORKDIR)
                 if os.path.isdir(os.path.join(WORKDIR, name))]

    maxBurstCount = max(distances, key=lambda d: len(list(glob("%s/%.3f/*" % (WORKDIR, d)))))
    runsToContinue = []
    for dist in distances:
        if len(list(glob("%s/%.3f/*" % (WORKDIR, dist)))) < maxBurstCount:
            runsToContinue.append(dist)
    if len(runsToContinue) > 0:
        log(GREEN + "Certain run sets require further burst runs. Starting total of %i runs.\n" %
            (len(runsToContinue) * args.burst) + END)
        for position in runsToContinue:
            newRunSet = RunSet(position)
            newRunSet.execute(IDs=newRunSet.restart())
        return False
    else:
        return True


def checkChiSquared():
    """Execute chi-squared analysis within a burst set"""
    log(BLUE + "Checking for chi-squared convergence of existing runs.\n" + END)

    # Get list of subdirectories
    distances = sorted([float(name.split('/')[-1]) for name in glob("%s/*.*" % WORKDIR)
                        if os.path.isdir(os.path.join(WORKDIR, name))], key=lambda n: float(n))

    runsToContinue = []
    errorMessage = RED + "Error: one or more dist_vs_t files are corrupt or missing.\n" + END

    for dist in distances:
        # Iterate over the existing runs.
        existingRuns = sorted([int(name) for name in os.listdir(WORKDIR + "/%.3f" % dist) if
                               os.path.isdir(os.path.join(WORKDIR + "/%.3f" % dist, name))])
        runsCount = len(existingRuns)
        chiSquaredSum = 0
        # Number of times already performed; each exclusion will leave out this many runs
        foldCount = int(round(runsCount / float(args.burst)))
        if runsCount % foldCount != 0:
            foldCount = 1
        log("Performing %i-fold calculations for distance " % foldCount + MAGENTA + "%.3f\n" % dist + END)

        def chunks(l, n):
            """Splits l into chunks of size n"""
            if n < 1:
                n = 1
            return [l[i:i + n] for i in range(0, len(l), n)]

        log("  ")
        for combination in chunks(existingRuns, foldCount):
            excludedInputList = []
            # Read data from all excluded runs
            for excludedBurstNum in combination:
                try:
                    with open(WORKDIR + "/%.3f" % dist + "/%i/dist_vs_t" % excludedBurstNum) as distVsT:
                        for line in distVsT:
                            excludedInputList.append(float(line.split()[1]))
                except IOError:
                    fail(errorMessage)
            # Find non-excluded runs
            checkDirectoriesList = sorted([int(name) for name in os.listdir(WORKDIR + "/%.3f" % dist) if
                                           os.path.isdir(os.path.join(WORKDIR + "/%.3f" % dist, name)) and int(name)
                                           not in combination])
            # Read included runs
            inputList = []
            for burstNum in checkDirectoriesList:
                try:
                    with open(WORKDIR + "/%.3f" % dist + "/%i/dist_vs_t" % burstNum) as distVsT:
                        for line in distVsT:
                            inputList.append(float(line.split()[1]))
                except IOError:
                    fail(errorMessage)

            def chiSquared(binListOne, binListTwo):
                """Compute the chi-squared distance between histograms. Returns the float distance."""

                # Compare starts of lists and add zeroes if necessary
                if binListOne[0][0] > binListTwo[0][0]:
                    addingAtIndex = 0
                    for thisTuple in binListOne:
                        if thisTuple[0] - binListTwo[addingAtIndex][0] > 0.009:
                            binListTwo.insert(addingAtIndex, (thisTuple[0], 0))
                            addingAtIndex += 1
                        else:
                            break
                elif binListTwo[0][0] > binListOne[0][0]:
                    addingAtIndex = 0
                    for thisTuple in binListTwo:
                        if thisTuple[0] - binListOne[addingAtIndex][0] > 0.009:
                            binListOne.insert(addingAtIndex, (thisTuple[0], 0))
                            addingAtIndex += 1
                        else:
                            break
                # Compare ends of lists and add zeroes if necessary
                if binListOne[-1][0] < binListTwo[-1][0]:
                    addingAtIndex = -1
                    for thisTuple in reversed(binListOne):
                        if binListTwo[addingAtIndex][0] - thisTuple[0] > 0.009:
                            if addingAtIndex == -1:
                                binListTwo.append((thisTuple[0], 0))
                            else:
                                binListTwo.insert(addingAtIndex + 1, (thisTuple[0], 0))
                            addingAtIndex -= 1
                        else:
                            break
                elif binListTwo[-1][0] < binListOne[-1][0]:
                    addingAtIndex = -1
                    for thisTuple in reversed(binListTwo):
                        if binListOne[addingAtIndex][0] - thisTuple[0] > 0.009:
                            if addingAtIndex == -1:
                                binListOne.append((thisTuple[0], 0))
                            else:
                                binListOne.insert(addingAtIndex + 1, (thisTuple[0], 0))
                            addingAtIndex -= 1
                        else:
                            break

                thisSum = 0
                for binNum in xrange(len(binListOne)):
                    a = binListOne[binNum][1]
                    b = binListTwo[binNum][1]
                    x = float((a - b) ** 2)
                    if a + b != 0:
                        thisSum += x / (a + b)
                return thisSum / 2

            def normalizedHist(dataList):
                """Wrapper for parallelization"""
                return histogram(dataList, binWidth=BINWIDTH, normalize=True)

            excludedHist, includedHist = parallelMap(normalizedHist, (excludedInputList, inputList))
            chiSquared = chiSquared(excludedHist, includedHist)
            if len(combination) == 1:
                log(str(combination[0]) + " ")
            else:
                log(str(combination) + " ")
            chiSquaredSum += chiSquared
        chiSquaredSum /= args.burst
        log("(average %.*f" % (args.precision, chiSquaredSum))
        if chiSquaredSum > args.chi_squared:
            # Continue this run
            log(", " + YELLOW + "exceeding" + END)
            runsToContinue.append(dist)
        log(")\n")
    if len(runsToContinue) > 0:
        log(GREEN + "Burst chi-squared analysis complete. Starting total of %i runs.\n" % (len(runsToContinue) *
                                                                                           args.burst) + END)
        for position in runsToContinue:
            newRunSet = RunSet(position)
            newRunSet.execute(IDs=newRunSet.restart())
        return False
    else:
        log(GREEN + "No restart runs necessary.\n" + END)
        return True


def analysis():
    """Execute the analysis process for burst runs involving histogram analysis"""
    global HISTPATH
    global STARTDIST
    global AVERAGEHILLWIDTH
    global DISTPATH
    log(CYAN + "Beginning analysis.\n" + END)

    # Get list of subdirectories
    distances = sorted([float(name.split('/')[-1]) for name in glob("%s/*.*" % WORKDIR)
                        if os.path.isdir(os.path.join(WORKDIR, name))])

    # Look for new runs
    log("Looking for new runs.\n")
    STARTDIST = float(min(glob("*.*"), key=lambda p: float(p)))
    dataList = []
    for thisFile in glob("*/*/dist_vs_t"):
        with open(thisFile) as infile:
            for line in infile:
                if line.split()[1] != "wham":
                    dataList.append(float(line.split()[1]))
    log("Loaded distance data.\n")

    # Generate histogram
    histogramData = list(reversed(histogram(dataList, binWidth=BINWIDTH)))
    histogramData = [item for item in histogramData if float(item[0]) < STARTDIST + MAXDIST + 10]  # Remove huge values
    log("Loaded histogram.\n")
    if args.hist is not None:
        with open(HISTPATH, 'w') as histFile:
            for b, c in histogramData:
                histFile.write("%.3f %i\n" % (b, c))

    # Begin traversing histogram output
    lastWasValley = True
    valleyEdges = []
    isFirstHill = True
    # Valley edges look like: [(last bin above threshold, first bin later above threshold),...]
    lastBin = ()
    firstHill = ()
    errorMessage = YELLOW + "Analysis ended: no bins fit the threshold. If runs should have been performed, the " \
                            "threshold parameter may be too big or too small.\n" + END

    for thisBin in histogramData:
        if lastBin == ():
            lastBin = thisBin
            continue

        # Check if it's a valley and if so append the first half of the tuple
        if thisBin[1] < args.thres <= lastBin[1]:
            if not lastWasValley:
                valleyEdges.append((lastBin,))
                lastWasValley = True

        # Check if it's a hill
        elif thisBin[1] >= args.thres > lastBin[1]:
            # Skip the first hill
            if isFirstHill:
                isFirstHill = False
                lastBin = thisBin
                lastWasValley = False
                firstHill = thisBin
                continue
            # Otherwise finish the tuple
            if lastWasValley:
                if len(valleyEdges[-1]) == 1:
                    valleyEdges[-1] = (valleyEdges[-1][0], thisBin)
                    lastWasValley = False
        lastBin = thisBin

    if len(valleyEdges) == 0 or len(valleyEdges[0]) != 2:
        fail(errorMessage)

    # Get the last hill
    lastHill = (valleyEdges[-2][1], valleyEdges[-1][0])

    # Invalidate valleys next to each other with no hill in between
    newValleyEdges = []
    skipNext = False
    for valleyNum in xrange(len(valleyEdges) - 1):
        if skipNext:
            skipNext = False
            continue
        thisValley = valleyEdges[valleyNum]
        nextValley = valleyEdges[valleyNum + 1]
        if thisValley[1][0] == nextValley[0][0]:
            # we have a conflict; keep the wider valley
            if thisValley[1][0] - thisValley[0][0] > nextValley[1][0] - nextValley[0][0]:
                # keep this one
                newValleyEdges.append(thisValley)
                skipNext = True  # Do not include the next one
                continue
            else:
                # do not keep this one; keep the other one
                continue
        else:  # no conflict
            newValleyEdges.append(thisValley)
    valleyEdges = newValleyEdges

    # Generate hills and valley widths
    hills = []
    valleyWidths = []
    for valley in valleyEdges:
        if len(valley) > 1:
            valleyWidths.append(valley[1][0] - valley[0][0])
        if len(hills) == 0:
            hills.append((firstHill,))
        hills[-1] = (hills[-1][0], valley[0])
        if len(valley) > 1:
            hills.append((valley[1],))

    # Strip anything
    if len(hills) == 0:
        fail(errorMessage)
    if len(hills[0]) != 2:
        hills.remove(hills[0])
    if len(hills) == 0:
        fail(errorMessage)
    if len(hills[-1]) != 2:
        hills.pop()
    if len(valleyEdges[-1]) != 2:
        valleyEdges.pop()

    hills.append(lastHill)

    # Generate average hill widths for the first time only
    sumHillWidths = 0
    if len(hills) < 2:
        fail(errorMessage)
    hillWidths = []
    for hill in hills:
        sumHillWidths += (hill[1][0] - hill[0][0])
        hillWidths.append(hill[1][0] - hill[0][0])
    if AVERAGEHILLWIDTH == 0:
        AVERAGEHILLWIDTH = sumHillWidths / len(hills)
    log("Average hill width is " + str(AVERAGEHILLWIDTH) + '\n')

    # Calculate number of simulations necessary in each valley
    simulationsPerValley = []
    for i in xrange(len(valleyWidths)):
        if i < len(hillWidths) and hillWidths[i] >= 0.001:  # Preferred
            simulationsPerValley.append(math.ceil(valleyWidths[i] / hillWidths[i]))
        else:  # Just use the average
            simulationsPerValley.append(math.ceil(valleyWidths[i] / AVERAGEHILLWIDTH))

    # Find original starting points to use instead of valleys to counter biasing
    newValleyEdges = []
    newValleyWidths = []
    if len(valleyEdges) + 1 == len(distances):
        for valleyNum in xrange(len(valleyEdges)):
            if valleyNum < len(hillWidths) - 1:
                newValleyLeft = distances[valleyNum] + 0.5 * hillWidths[valleyNum]
                newValleyRight = distances[valleyNum + 1] - 0.5 * hillWidths[valleyNum + 1]
            else:
                newValleyLeft = distances[valleyNum] + 0.5 * AVERAGEHILLWIDTH
                newValleyRight = distances[valleyNum] + 0.5 * AVERAGEHILLWIDTH + valleyWidths[valleyNum]
            newValleyEdges.append((newValleyLeft, newValleyRight))
            newValleyWidths.append(newValleyRight - newValleyLeft)
    else:
        for valleyNum in xrange(len(valleyEdges)):
            newValleyEdges.append((valleyEdges[valleyNum][0][0], valleyEdges[valleyNum][1][0]))
            newValleyWidths.append(valleyEdges[valleyNum][1][0] - valleyEdges[valleyNum][0][0])

    # Calculate starting positions.
    startPositions = []
    for valleyNum in xrange(len(valleyEdges)):
        valleyLeft = newValleyEdges[valleyNum][0]
        position = valleyLeft
        for newPositionNum in xrange(1, int(math.ceil(simulationsPerValley[valleyNum])) + 1):
            position += newValleyWidths[valleyNum] / (simulationsPerValley[valleyNum] + 1)
            if distances[0] < position < distances[-1]:
                startPositions.append(position)

    # Check if run is already in folder
    for position in startPositions:
        if "%.3f" % position in distances:
            startPositions.remove(position)

    log(BLUE + "The following runs are necessary: " + END)
    for position in startPositions:
        log(MAGENTA + "%.3f " % position + END)

    log(GREEN + "\nBurst analysis complete. Starting total of %i runs.\n" % len(startPositions) + END)
    if len(startPositions) > 0:
        with open(DISTPATH, 'r') as dist_vs_t:
            frames = dist_vs_t.read().splitlines()
        start = 0
        for position in startPositions:
            newRunSet = RunSet(position)
            runIDs, start = newRunSet.create(frames, start)
            newRunSet.execute(IDs=runIDs)
    elif args.normalize:
        log(BLUE + "Checking if run sets have the same number of burst runs.\n" + END)
        normalizePassed = checkNormalized()
        if normalizePassed:
            log(GREEN + "All runs completed, have equal runs, and meet threshold.\n" + END)
            sys.exit(0)
    else:
        log(GREEN + "All runs completed and meet threshold.\n" + END)
        sys.exit(0)


def checkFinishedRuns():
    """Check if all runs have finished execution."""
    finished = True

    def retry(path):
        log(YELLOW + UNDERLINE + "Warning:" + END + YELLOW + " run at %s has failed and will be restarted.\n" %
            path + END)
        with directory(path):
            system("rm PB* umb.out coord.nc end.rst mdinfo dist_vs_t finished executed > /dev/null 2> /dev/null")
            system("./run.sh > /dev/null 2> /dev/null")

    for runPath in [name for name in glob("*.*/*") if os.path.isdir(os.path.join(WORKDIR, name))]:  # Each run
        if glob(runPath + "/PB*o*"):
            output = open(glob(runPath + "/PB*o*")[0], 'r').read()
            if "unspecified launch failure" in output or "busy or unavailable" in output:
                # Try it again, since velocities will be regenerated
                retry(runPath)
                finished = False  # Not finished
                continue

        elif not os.path.isfile(runPath + "/finished"):
            if not os.path.isfile(runPath + "/executed"):  # It never even executed
                retry(runPath)
            finished = False  # Not finished
            continue

        elif (not os.path.isfile(runPath + "/coord.nc")) or os.stat(runPath + "/coord.nc")[6] <= 22 + len("coord.nc"):
            # Try it again
            retry(runPath)
            finished = False  # Not finished
            continue

    return finished


def eventLoop():
    """Handles pausing and analysis when runs complete."""
    global PAUSE
    global WORKDIR

    while True:
        # Check currently running processes
        finishedRuns = checkFinishedRuns()
        if finishedRuns:
            # Check chi-squared
            if args.chi_squared > 0:
                chiSquaredPassed = checkChiSquared()
                if chiSquaredPassed:
                    analysis()
            else:
                analysis()
        else:
            for i in xrange(PAUSE):
                time = (PAUSE - i)
                log((BOLD + "\rWill check again in " + MAGENTA + "%i" % time + END + BOLD +
                     " second" + ("s." if time > 1 else ".")).ljust(getTerminalWidth()) + END, toFile=False)
                sleep(1)
            log("\r")

    sys.exit(0)  # For safety; nothing should happen after the loop terminates.


if __name__ == "__main__":
    main()
