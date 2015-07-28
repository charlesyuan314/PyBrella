#!/bin/env python
# create_meta.py
# Generate the WHAM input files for PyBrella given a complete PyBrella run

from __future__ import division, print_function
from argparse import ArgumentParser, SUPPRESS
from glob import glob
from shared import *
from shutil import rmtree
from stat import S_IEXEC


parser = ArgumentParser(description="create WHAM input files for PyBrella")
parser.add_argument("--force", '-f', help="AMBER force constant (default is 5)", type=float, action="store", default=5)
parser.add_argument("--out", '-o', help="output directory", type=str, action=FullPath, default="wham")
parser.add_argument("--umbrella", '-u', help="PyBrella output directory", type=str, action=FullPath, default="umbrella")
parser.add_argument("--cont", help=SUPPRESS, action="store_true")
args = parser.parse_args()

WORKDIR = os.getcwd()
OUTPATH = args.out
if not os.path.isabs(args.out):
    OUTPATH = WORKDIR + "/" + args.out
UMBRELLAPATH = args.umbrella
if not os.path.isabs(args.umbrella):
    UMBRELLAPATH = WORKDIR + "/" + args.umbrella

if not args.cont:
    try:
        rmtree(OUTPATH)
    except OSError:
        pass
    os.mkdir(OUTPATH)

dists = [name.split('/')[-1] for name in glob("%s/*.*" % UMBRELLAPATH)
         if os.path.isdir(os.path.join(UMBRELLAPATH, name))]

dists.sort(key=lambda x: [int(y) for y in x.split('.')])
distFileNames = []
minPotentials = []
force = args.force * 2
sys.stderr.write("Processing distances: ")


def processDirectory(d):
    sys.stderr.write(d + " ")
    runs = [n.split('/')[-1] for n in glob("%s/%s/*" % (UMBRELLAPATH, d))
            if os.path.isdir(os.path.join(UMBRELLAPATH, name))]
    system("rm '%s/%s/dist_vs_t'" % (UMBRELLAPATH, d))
    for runDir in runs:
        system("cat '%s/%s/%s/dist_vs_t' >> '%s/%s/dist_vs_t'\n" % (UMBRELLAPATH, d, runDir, UMBRELLAPATH, d))
    system("mv %s/%s/dist_vs_t %s/dist_vs_t_%s" % (UMBRELLAPATH, d, OUTPATH, d))


if not args.cont:
    parMap(processDirectory, dists, n=(cpu_count() / 2))

sys.stderr.write("\nProcessing file lengths.\n")
maxLength = 0
distVsTs = glob("%s/dist_vs_t*" % OUTPATH)
for i, fp in enumerate(distVsTs):
    length = int(system("wc -l " + str(fp)).strip().split()[0])
    if length > maxLength:
        maxLength = length
    sys.stderr.write("\r" + str(int(float(i) / len(distVsTs) * 50)) + "% complete.")
for i, fp in enumerate(distVsTs):
    length = int(system("wc -l " + str(fp)).strip().split()[0])
    reps = int(maxLength / length)
    for _ in xrange(reps):
        system("cat %s >> %s_new" % (fp, fp))
    system("mv %s_new %s" % (fp, fp))
    sys.stderr.write("\r" + str(int(float(i) / len(distVsTs) * 50 + 50)) + "% complete.")
sys.stderr.write("\n")

sys.stderr.write("Creating meta.dat.\n")
for dist in dists:
    distFileNames.append("dist_vs_t_" + dist)
    distRst = open("%s/%s/0/dist.rst" % (UMBRELLAPATH, dist))
    for line in distRst:
        minPotentials.append(float(line.split()[4][3:-1]))
os.chdir(OUTPATH)
with open("meta.dat", 'w') as metaDat:
    for i in xrange(len(distFileNames)):
        metaDat.write("%s %.3f %.14f\n" % (distFileNames[i], minPotentials[i], force))

with open("run.sh", 'w') as runScript:
    runScript.write("wham %.3f %.3f 200 0.01 300 0 meta.dat result.dat\n" % (float(dists[0]), float(dists[-1])))

os.chmod("run.sh", os.stat("run.sh").st_mode | S_IEXEC)
os.chdir(WORKDIR)
