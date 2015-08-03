#!/bin/env python
# atom_tools.py
# Collection of PDB and atom/protein related functions

from numpy import array, dot, transpose, zeros, sum, sqrt
from numpy.linalg import svd, det
from shared import *

__author__ = 'Charles Yuan'


def atomDist(pdbLines, atom1, atom2):
    """Compute the distance for the PDB file between atom1 and atom2."""
    coord1 = []
    coord2 = []
    for pdbLine in pdbLines:
        # Get atoms
        if pdbLine.split()[0] == "ATOM":
            atomID = int(pdbLine.split()[1])
            if atomID == int(atom1):
                coord1 = array([float(pdbLine[(30 + i * 8):(38 + i * 8)]) for i in range(3)])
            elif atomID == int(atom2):
                coord2 = array([float(pdbLine[(30 + i * 8):(38 + i * 8)]) for i in range(3)])
    diff = coord1 - coord2
    return sqrt(dot(diff, diff))  # vector distance


def calcCenter(pdbLines):
    """Calculates weighted center of structure (center of mass) based on atomic masses.
    Originally part of pdbTools by Mike Harms."""

    # Weights for atoms
    ATOM_WEIGHTS = {"H": 1.00794,
                    "D": 2.01410178,  # deuterium
                    "HE": 4.00,
                    "LI": 6.941,
                    "BE": 9.01,
                    "B": 10.811,
                    "C": 12.0107,
                    "N": 14.0067,
                    "O": 15.9994,
                    "F": 18.998403,
                    "NE": 20.18,
                    "NA": 22.989769,
                    "MG": 24.305,
                    "AL": 26.98,
                    "SI": 28.09,
                    "P": 30.973762,
                    "S": 32.065,
                    "CL": 35.453,
                    "AR": 39.95,
                    "K": 39.0983,
                    "CA": 40.078,
                    "SC": 44.96,
                    "TI": 47.87,
                    "V": 50.94,
                    "CR": 51.9961,
                    "MN": 54.938045,
                    "FE": 55.845,
                    "CO": 58.93,
                    "NI": 58.6934,
                    "CU": 63.546,
                    "ZN": 65.409,
                    "GA": 69.72,
                    "GE": 72.64,
                    "AS": 74.9216,
                    "SE": 78.96,
                    "BR": 79.90,
                    "KR": 83.80,
                    "RB": 85.47,
                    "SR": 87.62,
                    "Y": 88.91,
                    "ZR": 91.22,
                    "NB": 92.91,
                    "MO": 95.94,
                    "TC": 98.0,
                    "RU": 101.07,
                    "RH": 102.91,
                    "PD": 106.42,
                    "AG": 107.8682,
                    "CD": 112.411,
                    "IN": 114.82,
                    "SN": 118.71,
                    "SB": 121.76,
                    "TE": 127.60,
                    "I": 126.90447,
                    "XE": 131.29,
                    "CS": 132.91,
                    "BA": 137.33,
                    "PR": 140.91,
                    "EU": 151.96,
                    "GD": 157.25,
                    "TB": 158.93,
                    "W": 183.84,
                    "IR": 192.22,
                    "PT": 195.084,
                    "AU": 196.96657,
                    "HG": 200.59,
                    "PB": 207.2,
                    "U": 238.03}

    # Calculate the center of mass of the protein (assuming all atoms of one type have the
    # same mass).
    coord = []
    masses = []
    for l in pdbLines:
        # Skip non ATOM lines
        if l[0:4] != "ATOM":
            continue

        # Grab coordinates
        coord.append([float(l[(30 + i * 8):(38 + i * 8)]) for i in range(3)])

        # Grab mass of each atom.  First try to grab atom type entry.
        # If it's missing, guess from the atom name column.
        atom_type = l[73:].strip()
        if atom_type == "":
            if l[12] == " ":
                atom_type = l[13]
            elif l[12] == "H":
                atom_type = "H"
            else:
                atom_type = l[12:14].strip()

        # Append the atom mass
        if atom_type == "C" or atom_type not in ATOM_WEIGHTS:
            masses.append(12.0107)
        else:
            masses.append(ATOM_WEIGHTS[atom_type])

    num_atoms = len(coord)

    total_mass = sum(masses)
    weights = [m / total_mass for m in masses]
    center = array([sum([coord[i][j] * weights[i] for i in range(num_atoms)]) for j in range(3)])

    return center


def calcCenterAtoms(pdbLines):
    """Given a pdb file, calculate the centers of the ligand and the protein"""
    pPDB = []
    lPDB = []
    for line in pdbLines:
        if line[0:4] != "ATOM":
            continue
        if line[17:20] in RESIDUES:
            pPDB.append(line)
        elif line[17:20] != "WAT":
            lPDB.append(line)
    pCenter = calcCenter(pPDB)
    lCenter = calcCenter(lPDB)

    lAtom, lAtomCoord = closestAtom(lPDB, lCenter)
    pAtom, pAtomCoord = closestAtom(pPDB, pCenter)
    return lAtom, pAtom, lAtomCoord, pAtomCoord, lCenter, pCenter


def closestAtom(pdbLines, coords):
    """Given a NumPy array of three-dimensional coordinates,
    return (ID, coord) of the closest atom in the file pdb."""

    # Calculate the coordinates of each atom
    coord = []  # {(atomID, coord), ...}
    for l in pdbLines:
        # Skip non ATOM lines
        if l[0:4] != "ATOM":
            continue
        # Grab coordinates
        coord.append((int(l.split()[1]), array([float(l[(30 + i * 8):(38 + i * 8)]) for i in range(3)])))

    def c(atom):
        diff = atom[1] - coords
        return sqrt(dot(diff, diff))
    sortedAtoms = min(coord, key=c)

    return sortedAtoms


def rmsdDist(pdbLines, refCoords, segments=None):
    """Compute the distance for the PDB file and the reference file."""

    def rmsd(P, Q):
        """Returns RMSD between 2 sets of nx3 numpy array using the Kabsch algorithm.
           Modified from rmsd by Jimmy Charnley Kromann.
           License at https://github.com/charnley/rmsd/blob/master/LICENSE"""
        def kabsch_rmsd(x, y):
            C = dot(transpose(x), y)
            V, S, W = svd(C)
            if (det(V) * det(W)) < 0.0:
                S[-1] = -S[-1]
                V[:, -1] = -V[:, -1]
            U = dot(V, W)
            x = dot(x, U)
            D = len(x[0])
            N = len(x)
            r = 0.0
            for v, w in zip(x, y):
                r += sum([(v[k] - w[k]) ** 2.0 for k in range(D)])
            return sqrt(r / N)
        # Calculate RMSD for centered sets; just subtract the centroids for P and Q
        # noinspection PyTypeChecker
        return kabsch_rmsd(P - sum(P) / len(P), Q - sum(Q) / len(Q))

    processedLines = []
    if segments is None:  # Process everything
        processedLines = [pdbLine for pdbLine in pdbLines if pdbLine[13:15] == "CA"]
    else:  # Process only segments
        for pdbLine in pdbLines:
            if pdbLine[13:15] == "CA":
                resID = int(pdbLine[21:28])
                for tup in segments:
                    if tup[0] <= resID <= tup[1]:
                        processedLines.append(pdbLine)

    coords = zeros((len(processedLines), 3), float)
    for index, line in enumerate(processedLines):
        for j in range(3):
            coords[index, j] = float(line[(30 + j * 8):(38 + j * 8)])

    return rmsd(coords, refCoords)
