#****************************************************************************
# Copyright (C) 2013 SUNCAT
# This file is distributed under the terms of the
# GNU General Public License. See the file `COPYING'
# in the root directory of the present distribution,
# or http://www.gnu.org/copyleft/gpl.txt .
#****************************************************************************


import numpy as np
from ase import constraints


class specobj:
    #small species class containing the attributes of a species
    def __init__(self, s='X', mass=0., magmom=0., U=0., J=0., U_alpha=0.):
        self.s = s
        self.mass = mass
        self.magmom = magmom
        self.U = U
        self.J = J
        self.U_alpha = U_alpha


#add 'd0' to floating point number to avoid random trailing digits in Fortran input routines
def num2str(x):
    s = str(x)
    if s.find('e')<0:
        s += 'd0'
    return s

#convert python to fortran logical
def bool2str(x):
    if x:
        return '.true.'
    else:
        return '.false.'

#convert some of ase's constraints to pw.x constraints for pw.x internal relaxation
#returns constraints which are simply expressed as setting force components
#as first list and other contraints that are implemented in espresso as
#second list
def convert_constraints(atoms):
    if atoms.constraints:
        n = len(atoms)
        if n==0:
            return [],[]
        forcefilter = []
        otherconstr = []
        for c in atoms.constraints:
            if isinstance(c, constraints.FixAtoms):
                if len(forcefilter)==0:
                    forcefilter = np.ones((n,3), np.int)
                forcefilter[c.index] = [0,0,0]
            elif isinstance(c, constraints.FixCartesian):
                if len(forcefilter)==0:
                    forcefilter = np.ones((n,3), np.int)
                forcefilter[c.a] = c.mask
            elif isinstance(c, constraints.FixBondLengths):
                for d in c.constraints:
                    otherconstr.append("'distance' %d %d" % (d.indices[0]+1,d.indices[1]+1))
            elif isinstance(c, constraints.FixBondLength):
                otherconstr.append("'distance' %d %d" % (c.indices[0]+1,c.indices[1]+1))
            elif isinstance(c, constraints.FixInternals):
            # we ignore the epsilon in FixInternals because there can only be one global
            # epsilon be defined in espresso for all constraints
                for d in c.constraints:
                    if isinstance(d, constraints.FixInternals.FixBondLengthAlt):
                        otherconstr.append("'distance' %d %d %s" % (d.indices[0]+1,d.indices[1]+1,num2str(d.bond)))
                    elif isinstance(d, constraints.FixInternals.FixAngle):
                        otherconstr.append("'planar_angle' %d %d %d %s" % (d.indices[0]+1,d.indices[1]+1,d.indices[2]+1,num2str(np.arccos(d.angle)*180./np.pi)))
                    elif isinstance(d, constraints.FixInternals.FixDihedral):
                        otherconstr.append("'torsional_angle' %d %d %d %d %s" % (d.indices[0]+1,d.indices[1]+1,d.indices[2]+1,d.indices[3]+1,num2str(np.arccos(d.angle)*180./np.pi)))
                    else:
                        raise NotImplementedError, 'constraint '+d.__name__+' from FixInternals not implemented\n' \
                            'consider ase-based relaxation with this constraint instead'                        
            else:
                raise NotImplementedError, 'constraint '+c.__name__+' not implemented\n' \
                    'consider ase-based relaxation with this constraint instead'
        return forcefilter,otherconstr
    else:
        return [],[]
