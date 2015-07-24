#!/bin/env python
# shared.py
# Collection of utilities for PyBrella, lineMD, and related scripts

from __future__ import division, print_function
from argparse import Action
from contextlib import contextmanager
from itertools import izip
import math
from multiprocessing import cpu_count, Pipe, Process, Queue
from operator import itemgetter
import os
from subprocess import PIPE, Popen, STDOUT
import sys


__author__ = 'Charles'


"""Define ANSI sequences"""
BOLD = "\033[1m"
UNDERLINE = "\033[4m"
BLACK = "\033[30m"
RED = "\033[31m"
GREEN = "\033[32m"
YELLOW = "\033[33m"
BLUE = "\033[34m"
MAGENTA = "\033[35m"
CYAN = "\033[36m"
WHITE = "\033[37m"
END = "\033[0m"


RESIDUES = ["ACE", "NHE", "NME", "ALA", "ARG", "ASN", "ASP", "CYS", "CYX", "GLU",
            "GLN", "GLY", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "MET", "PHE",
            "PRO", "SER", "THR", "TRP", "TYR", "VAL"]


@contextmanager
def directory(d):
    """Context manager that executes a block of code inside directory d and returns to the
    previous working directory when done."""
    curdir = os.getcwd()
    try:
        os.chdir(d)
        yield
    finally:
        os.chdir(curdir)


class FullPath(Action):
    """Expand user and relative paths."""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


def cmp_to_key(mycmp):
    """Convert a cmp= function into a key= function"""

    class K(object):
        def __init__(self, obj):
            self.obj = obj

        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0

        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0

        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0

        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0

        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0

        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0

    return K


def getTerminalWidth():
    """Code snippet to find terminal width in characters."""
    env = os.environ

    def ioctl_GWINSZ(f):
        try:
            from fcntl import ioctl
            from struct import unpack
            from termios import TIOCGWINSZ

            c = unpack('hh', ioctl(f, TIOCGWINSZ, "1234"))
        except StandardError:
            return
        return c

    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            # noinspection PyTypeChecker
            os.close(fd)
        except StandardError:
            pass
    if not cr:
        cr = (env.get("LINES", 25), env.get("COLUMNS", 80))
    return int(cr[1])


def fail(string):
    """Print an error and exit."""
    log('\n' + string)
    sys.exit(1)


def system(command):
    """Execute a command on the shell. Blocks the thread and returns output. Do not use for commands
    based on user input."""
    return Popen(command, stderr=STDOUT, stdout=PIPE, shell=True).communicate()[0]


def unblocked_system(command):
    """Execute a command on the shell. Does not block the thread or return output.
    Do not use for commands based on user input."""
    Popen(command, shell=True)


def compress(name):
    """Compress a file with gzip."""
    system("gzip -f " + str(name))


def decompress(name):
    """Decompress a file with gzip."""
    system("gunzip -f " + str(name))


def parallelMap(f, x):
    """Execute function f for each element in set x. Executes functions defined
    anywhere and parallelizes with multiprocessing. Returns a list of outputs
    and blocks the main process until complete."""
    # Warning: will spawn len(X) processes. Use parMap or splitMap for large X
    pipe = []
    for _ in x:
        pipe.append(Pipe())

    def spawn(f2):
        def fun(p2, x2):
            p2.send(f2(x2))
            p2.close()

        return fun

    proc = [Process(target=spawn(f), args=(c, k)) for k, (p, c) in izip(x, pipe)]
    [p.start() for p in proc]
    [p.join() for p in proc]
    return [p.recv() for (p, c) in pipe]


def parMap(f, X, n=(cpu_count() / 2), silent=False):
    """Executes f over X with nprocs processes at once. Superior to
    splitMap because items close in X execute close to each other in time."""
    # splitMap splits the list, resulting in distant tasks running at once.

    n = int(n)

    def chunks(seq, k):
        """Returns n subsequences of roughly equal length from seq, useful for parmap.
        seq is divided into k lists not by splitting but by item-by-item distribution."""
        results = []
        for _ in range(k):
            results.append([])
        for j, e in enumerate(seq):
            results[j % k].append(e)
        return results

    testLists = list(chunks(X, n))

    def m(func, l):
        """Flat sequential map function that writes status"""
        o = []
        tot = len(l)
        for j, t in enumerate(l):
            if not silent:
                log(("\rApproximately %.3f%% complete" % (100 * float(j) / tot))
                    .ljust(getTerminalWidth()))
            o.append((t, func(t)))
        if not silent:
            log("\r")
        return o

    q_in = Queue(1)
    q_out = Queue()

    def fun(f_in, m_func, qin, qout):
        """Wrapper function for parallelization"""
        while True:
            j, x_in = qin.get()
            if j is None:
                break
            qout.put((j, m_func(f_in, x_in)))

    proc = [Process(target=fun, args=(f, m, q_in, q_out)) for _ in range(n)]
    for p in proc:
        p.daemon = True
        p.start()

    sent = [q_in.put((i, x)) for i, x in enumerate(testLists)]
    [q_in.put((None, None)) for _ in range(n)]
    res = [q_out.get() for _ in range(len(sent))]

    [p.join() for p in proc]
    output = sorted([item for sublist in [x for i, x in sorted(res)] for item in sublist],
                    key=itemgetter(0))
    return [item[1] for item in output]


def splitMap(f, X, n=(cpu_count() / 2)):
    """Execute function f for each element in set X, running n operations at once. Executes
    functions defined anywhere and parallelizes with multiprocessing. Returns a list of
    outputs and blocks the main process until complete."""

    n = int(n)

    def chunks(l, k):
        """Yield successive k-sized chunks from l."""
        for i in xrange(0, len(l), k):
            yield l[i:i + k]

    output = []
    for chunk in list(chunks(X, n)):
        output += parallelMap(f, chunk)
    return output


def log(message, toFile=False):
    """Log a message. toFile is ignored."""
    if toFile:
        pass
    sys.stderr.write(message)


def frange(limit1, limit2=None, increment=1.):
    """
    Range function that accepts floats (and integers).
    The returned value is an iterator.  Use list(frange) for a list.
    """

    if limit2 is None:
        limit2, limit1 = limit1, 0.
    else:
        limit1 = float(limit1)

    count = int(math.ceil(limit2 - limit1) / increment)
    return (limit1 + n * increment for n in range(count))
