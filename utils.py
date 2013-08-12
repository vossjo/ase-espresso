import numpy as np
from ase import constraints

def uniqueness(list1, list2):
    """
    Find the elements that are belonging to a group in both list1 and list2,
    where a group is defines as 2 or more elements that share the same
    value. 
    """
    assert len(list1) == len(list2)
    if len(list1)==1:
        return [1]
    l1_u = np.unique(list1)
    l2_u = np.unique(list2)
    unique = np.zeros(len(list1))
    kk = 0
    for u_l1 in l1_u:
        UK1 = np.where(list1 == u_l1)[0]
        UL1 = [pp in UK1 for pp in range(len(list1))]
        for u_l2 in l2_u:
            UK2 = np.where(list2 == u_l2)[0]
            UL2 = [pp in UK2 for pp in range(len(list1))]
            UUL = [UL1[pp]*UL2[pp] for pp in range(len(list1))]
            if len(np.where(np.array(UUL) != 0)[0]) == 0:
                continue 
            kk += 1
            unique += [kk*UUL[pp] for pp in range(len(list1))]
    # fill out zeros
    umax = np.max(unique)
    zeros = np.where(unique==0)[0]
    for ppk, pp in enumerate(zeros):
        unique[pp] += ppk+umax
    return unique.astype(int)

#add 'd0' to floating point number to avoid random trailing digits in Fortran input routines
def num2str(x):
    s = str(x)
    if s.find('e')<0:
        s += 'd0'
    return s


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
