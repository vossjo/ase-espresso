#!/usr/local/bin/python


#****************************************************************************
# Copyright (C) 2013 SUNCAT
# This file is distributed under the terms of the
# GNU General Public License. See the file `COPYING'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#****************************************************************************


from sys import argv, exit, stderr
if len(argv) != 3:
    print >>stderr, 'usage: ' + argv[0] + ' pw-log output.traj'
    exit(1)

import numpy as np
from ase.io import Trajectory
from ase import Atoms
from ase.calculators.singlepoint import SinglePointDFTCalculator
from os import popen

import ase.units
hartree = ase.units.Hartree
rydberg = ase.units.Rydberg
bohr = ase.units.Bohr

# Get energy from log.
def get_total_energy(s):
    a = s.readline()
    while a != '' and a[:17] != '!    total energy':
        a = s.readline()
    if a == '':
        return None
    energy = float(a.split()[-2]) * rydberg
    for i in range(11):
        a = s.readline()
        if '     smearing contrib.' in a:
            break
    # correct for finite temperature entropy term
    # in case of finite temp. smearing
    if a[:22] == '     smearing contrib.':
        energy -= 0.5 * float(a.split()[-2]) * rydberg
    return energy

p = popen('grep -n Giannozzi ' + argv[1] + ' 2>/dev/null | tail -1', 'r')
try:
    n = int(p.readline().split()[0].strip(':'))
except:
    print >>stderr, 'No valid pw-log at ' + argv[1] + ' found.'
    p.close()
    exit(2)
p.close()

s = open(argv[1], 'r')
# skip over previous runs in log in case the current log has been
# appended to old ones
for i in range(n):
    s.readline()

a = s.readline()
while a[:11] != '     celldm':
    a = s.readline()
alat = float(a.split()[1]) / 1.889726
a = s.readline()
while a[:12] != '     crystal':
    a = s.readline()
cell = []
for i in range(3):
    cell.append([float(x) for x in s.readline().split()[3:6]])
cell = np.array(cell)
a = s.readline()
while a[:12] != '     site n.':
    a = s.readline()
pos = []
syms = ''
y = s.readline().split()
while len(y) > 0:
    nf = len(y)
    pos.append([float(x) for x in y[nf - 4:nf - 1]])
    syms += y[1].strip('0123456789')
    y = s.readline().split()
pos = np.array(pos) * alat
natoms = len(pos)

# create atoms object with coordinates and unit cell
# as specified in the initial ionic step in log
atoms = Atoms(syms, pos, cell=cell * alat, pbc=(1, 1, 1))

# Attach calculator object to atoms.
calc = SinglePointDFTCalculator(atoms)
atoms.set_calculator(calc)
# Get total energy at first ionic step and store it in the calculator.
en = get_total_energy(s)
if en is not None:
    calc.results['energy'] = en
else:
    print >>stderr, 'no total energy found'
    exit(3)

# Write to .traj file
traj = Trajectory(argv[2], 'w')
traj.write(atoms)

# Append the following images.
a = s.readline()
while a != '':
    while a[:7] != 'CELL_PA' and a[:7] != 'ATOMIC_' and a != '':
        a = s.readline()
    if a == '':
        break
    if a[0] == 'A':
        coord = a.split('(')[-1]
        for i in range(natoms):
            pos[i][:] = s.readline().split()[1:4]

        if coord == 'alat)':
            atoms.set_positions(pos * alat)
        elif coord == 'bohr)':
            atoms.set_positions(pos * bohr)
        elif coord == 'angstrom)':
            atoms.set_positions(pos)
        else:
            atoms.set_scaled_positions(pos)
        # Get total energy at 2nd,3rd,...nth ionic step.
        en = get_total_energy(s)
        # Initialize calculator with updated atoms.
        calc = SinglePointDFTCalculator(atoms)
        atoms.set_calculator(calc)
        calc.results['energy'] = en

        if en is not None:
            traj.write(atoms)
        else:
            break
    else:
        for i in range(3):
            cell[i][:] = s.readline().split()
        atoms.set_cell(cell * alat, scale_atoms=False)
    a = s.readline()
