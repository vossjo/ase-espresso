"""
This module contains functionality for the QUANTUM ESPRESSO calculator
in order to:

   - read the results in text and xml format
   - read a .in file and return all the data
   - build an Atoms objects starting from an .in file


routines:


def read_quantumespresso_textoutput(filename):
def read_quantumespresso_xmloutput(filename):

def write_atoms(infilename, myprefix="dry_pwscf", verbose = False):
def read_in_file(infilename):


"""

import os
import string
import copy

import sys
try:
    import xml.etree.cElementTree as ET
except ImportError:
    try:
        import xml.etree.ElementTree as ET
    except ImportError:
        print("unable to import ElementTree")

import numpy as np
from ase.atoms import Atoms
import ase.units
hartree = ase.units.Hartree
rydberg = ase.units.Rydberg
bohr = ase.units.Bohr

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
        correction = 0.5 * float(a.split()[-2]) * rydberg
        energy -= correction
    return energy


# dummy calculator
class calculator:

    def set_energy(self, e):
        self.energy = e

    def get_potential_energy(self):
        return self.energy

    def get_calculator(self):
        return self

    def notimpl(self, apply_constraint=False):
        raise NotImplementedError


# -------------------------------------------------------------------------------------------------------------------------------
#
#  read text output
#
# -------------------------------------------------------------------------------------------------------------------------------    
        
def read_quantumespresso_textoutput(filename, verbose=False):
    """
    read the output text file from a QE run and extracts physical quantities
    returns a dictionary with the quantities
    """
    
    alltext = open(filename).read().lower()                             # reads all the lines in the text file and converts them in
                                                                        # lowercase
    
    ## [RFC] some 'assert' statements about NOT having some error string in the text
    
    alllines = alltext.split('\n')                                      # creates a list of lines from the text stream

        #kpts_num = None
        #nbands = None
        #nelect = None
    _no_spin_orbit = None
    Results = {}
    

    Results['error'] = False    
    Results['error_message'] = []
    Results['ExchangeCorrelation'] = None

    Results['error_found'] = False

    matches=[0]                                                         #  the list will store the line numbers at which each iteration ends
                                                                        #  in that way, in case of multiple iterations, we will be allowed to
                                                                        #  excerpt the last one
    for n, l in enumerate(alllines):                                    #  [RFC] are there other useful infos in the preambole?
        
        if l.rfind("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%") > -1:
            if Results['error'] :
                Results['error_message'].append(l)
                break;
            Results['error'] = True
            
        if Results['error']:
            Results['error_message'].append(l)
                    
        if (l.strip()=='writing output data file pwscf.save'):
            matches.append(n)

        if l.rfind('number of atoms/cell') > -1:
            t = l.split()
            pos = t.index('=')
            Results['nat'] = int(t[pos+1])
            
        if l.rfind('number of electrons') > -1:
            t = l.split()
            pos = t.index('=')
            Results['nelec'] = float(t[pos+1])

        if l.rfind('kohn-sham states') > -1:                        # the number of bands
            t = l.split()
            pos = t.index('states=')
            Results['nbands'] = int(t[pos+1])

        if l.rfind('number of k points=') > -1:                     # the number of kpoints
            t = l.split()
            pos = t.index('points=')
            Results['kpts_num'] = int(t[pos+1])

        if l.rfind('Eexchange-correlation') > -1:
            t = l.split()
            pos = t.index('=')
            Results['exchange_correlation'] = t[pos+1:len(t)-1]

        if l.rfind('without spin-orbit') > -1:
            _no_spin_orbit = True
            
    if(len(matches)>1):
        last = len(matches)
        lines = alllines[matches[last-2]:matches[last-1]]
    else:
        lines = alllines
    
        
    alltext  = [];                                                      # some memory garbage collection would be good here
    alllines = [];
            
    Results['EnergyFound'] = False
    CAUGHT = -1

    if (Results['error']):
        return False

            
#        self.etotal = None;
#        self.etotal_accuracy = None;
    Results['niter'] = None
    Results['pressure'] = None
    Results['stress'] = None
    Results['total_magnetization'] = None
    Results['absolute_magnetization'] = None
    Results['fermi_energy'] = None
    Results['kpts'] = None
    Results['atoms_forces'] = None
#                    _old_cell = self.cell.copy()
    Results['cell'] = None
    Results['atomic_positions'] = None

    for n, l in enumerate(lines):                                       # loop on the lines from output

                                                                        # some quantities, like the stress tensor or the kpoints
                                                                        #  span several lines. we use the trick of tagging the
                                                                        #  begin of those lines by setting the variable "CAUGHT"
                                                                        #  A non-negative value means that we are still inside
                                                                        #  such a sub-block. The number of lines still left in the
                                                                        #  sub-block is given by the value of CAUGHT, which is
                                                                        #  decreased by one at each line.

            
        if CAUGHT < 0:
            # we are not in a sub-block
                                                                        # find Total Energy and accuracy
            if len(l)>0 and (l[0] == '!'):  # and l[1] == 'total' and l[2] == 'energy'):
                t = l.split();
                if(t[1] == 'total' and t[2] == 'energy'):
                    #print "ENERGY FOUND"
                    pos = t.index('ry')
                    Results['etotal'] = float(t[pos-1])                            # total energy in Rydberg
                    Results['EnergyFound'] = True

            elif Results['EnergyFound'] and l.rfind('scf accuracy') > -1:
                #print "ACCURACY"
                t = l.split()
                pos = t.index('ry')
                Results['etotal_accuracy'] = float(t[pos-1])

            if l.rfind('fermi energy') > -1:
                t = l.split()
                pos = t.index('is')
                Results['fermi_energy'] = float(t[pos+1])

            if l.rfind('magnetization') > -1:
                t = l.split()
                pos = t.index('magnetization')
                if(t[pos-1] == 'total'):
                    Results['total_magnetization'] = float(t[pos+2])
                elif(t[pos-1] == 'absolute'):
                    Results['absolute_magnetization'] = float(t[pos+2])
                    
            if l.rfind('convergence has been achieved in') > -1:
                t = l.split()
                pos = t.index('iterations')
                Results['niter'] = int(t[pos-1])

            if l.rfind('lattice parameter (alat)') > -1:
                t = l.split()
                pos = t.index('=')
                Results['alat'] = float(t[pos+1])
                
            if l.rfind('stress') > -1:                                  # find the stress tensor
                t = l.split()
                pos = t.index('stress')
                if (pos > 0 and t[pos-1]=='total'):
                    if l.rfind('p=') > -1:                              # find the pressure
                        pos = t.index('p=')
                        Results['pressure'] = float(t[pos+1])
                                                                        # set-up to manage the sub-block
                    CAUGHT = 3
                    CAUGHT_DISCARD=0
                    CAUGHT_WHAT = 'stress'
                    CAUGHT_LINES = []

            if l.rfind('number of k points=') > -1:                           # find the kpoints
                CAUGHT = Results['kpts_num']
                CAUGHT_DISCARD = 1                                      # discard 1 line at the beginning
                CAUGHT_WHAT = 'kpts'
                CAUGHT_LINES = []
                        
                        
            if l.rfind('new k-points:') > -1:                           # find the kpoints
                CAUGHT = Results['kpts_num']
                CAUGHT_DISCARD = 0                                      # discard 0 line at the beginning
                CAUGHT_WHAT = 'kpts'
                CAUGHT_LINES = []
            
            if l.rfind('forces acting on atoms') > -1:                  # find forces acting on atoms
                CAUGHT = Results['nat']
                CAUGHT_DISCARD = 1
                CAUGHT_WHAT = 'forces'
                CAUGHT_LINES = []

            if l.rfind('crystal axes:') > -1:
                CAUGHT = 3
                CAUGHT_DISCARD = 0
                CAUGHT_WHAT = 'old cell'
                CAUGHT_LINES = []                           
                
            if l.rfind('final estimate of lattice vectors') > -1:
                CAUGHT = 3
                CAUGHT_DISCARD = 0
                CAUGHT_WHAT = 'new cell'
                CAUGHT_LINES = []

            if l.rfind('atomic_positions') > -1:
                CAUGHT= Results['nat']
                CAUGHT_DISCARD = 0
                CAUGHT_WHAT = 'positions'
                CAUGHT_LINES = []
                
        elif CAUGHT > 0:                                                # We are inside a sub-blok
            if(CAUGHT_DISCARD == 0):                                    # gets the line but if this line should be
                CAUGHT_LINES.append(l)                                  # discarded instead
                CAUGHT -= 1
            else:
                CAUGHT_DISCARD -= 1

        else:                                                           # process the caught lines
                
            if(CAUGHT_WHAT == 'stress'):
                                                                        # get the stress tensor in kbar
                Results['stress']=[]
                number = len(CAUGHT_LINES)                              # i.e. 3, in this case
                for n in range(number):
                    Results['stress'].append([float(el) for el in CAUGHT_LINES[n].split()[3:6]])   #NOTE it is in kbar        
                    
            elif (CAUGHT_WHAT == 'kpts'):                               # get the kpoints
                Results['kpts']    = []
                Results['kpts_wk'] = []
                number = len(CAUGHT_LINES)
                for n in range(number):
                    try:
                        Results['kpts'].append([float(get_number(el)) for el in CAUGHT_LINES[n].split()[4:7]])
                        Results['kpts_wk'].append([float(CAUGHT_LINES[n].split()[9])])
                    except:
                        print("Trouble with " + CAUGHT_LINES[n])

            elif (CAUGHT_WHAT == 'forces'):
                Results['atoms_forces'] = []
                number = len(CAUGHT_LINES)
                for n in range(number):
                    Results['atoms_forces'].append([float(get_number(el)) for el in CAUGHT_LINES[n].split()[6:9]])

                # elif (CAUGHT_WHAT == 'old cell'):
                #     old_cell = []
                #     number = len(CAUGHT_LINES)
                #     for n in range(number):
                #         old_cell.append([float(el) for el in CAUGHT_LINES[n].split()[3:5]])

            elif (CAUGHT_WHAT == 'new cell'):
                Results['cell'] = []
                number = len(CAUGHT_LINES)
                for n in range(number):
                    Results['cell'].append([float(el) for el in CAUGHT_LINES[n].split()])
            
            elif (CAUGHT_WHAT == 'positions'):
                Results['atomic_positions'] = []
                number = len(CAUGHT_LINES)
                for n in range(number):
                    Results['atomic_positions'].append([float(el) for el in CAUGHT_LINES[n].split()[1:4]])        
                                        
            CAUGHT_LINES = None
            CAUGHT_WHAT = None
            CAUGHT = -1
                

                        
    if verbose:
        print "etotal is           : %s " % str(Results['etotal'])
        print "etotal accuracy is  : %s " % str(Results['etotal_accuracy'])
        print "niter is            : %s " % str(Results['niter'])
        print "nelect is           : %s " % str(Results['nelect'])
        print "Results['nbnds is           : %s " % str(Results['nbands'])
        print "Fermi energy is     : %s " % str(Results['fermi_energy'])
        print "tot magnetizations  : %s " % str(Results['total_magnetization'])
        print "abs magnetizations  : %s " % str(Results['absolute_magnetization'])
        print "ExchangeCorrelation : "
        if (not Results['ExchangeCorrelation'] is None):
            print "%s " % [str(el) for el in Results['ExchangeCorrelation']]
        else:
            print "None"
        if(not Results['pressure'] is None):
            print "Pressure is         : %s " % str(Results['pressure'])
        if(not Results['stress'] is None):  
            print "Stress is           : %s " % str(Results['stress'])
        print "kpoints found       : %s " % str(Results['kpts_num'])
        print "kpoints are         : "
        if(not Results['kpts'] is None):
            for n in range(len(Results['kpts'])):
                print "\t %s wk: %s" % (str(Results['kpts'][n]), str(Results['kpts_wk'][n]))
        if(not Results['atoms_forces'] is None):
            print "forces are          : "
            print Results['atoms_forces']
        print "alat is             : %s" % str(Results['alat'])
        if (not Results['cell'] is None):
            print "cell is         : "
            print Results['new_cell']
        if( not Results['atomic_positions'] is None):
            print "atomic positions    : "
            print Results['new_atomic_positions']

    return Results
                

def get_bands_for_kpoints(mylines, nk, nb):
    energy_array=[]
    for n, l in enumerate(mylines):

        if l.expandtabs(1).isspace():
            pass
        t = l.split()
        if(t[0] == 'k'):
            pass
        energy_array.append([float(el) for el in t])
        
        
    
def get_number(string):
    for t in string.split():
        array = []
        for c in t:
            if (c.isdigit() or (c in '.e+-')):
                array.append(str(c))
    return ''.join([c for c in array])
    





# -------------------------------------------------------------------------------------------------------------------------------
#
#  read xml output
#
# -------------------------------------------------------------------------------------------------------------------------------    

bravais_lattice_to_ibrav = {"free":0,
                            "cubic p (sc)": 1,
                            "cubic f (fcc)": 2,
                            "cubic i (bcc)": 3,
                            "hexagonal and trigonal p": 4,
                            "trigonal r": 5,
                            "tetragonal p (st)": 6,
                            "tetragonal i (bct)": 7,
                            "orthorhombic p": 8,
                            "orthorhombic base-centered(bco)": 9,
                            "orthorhombic face-centered": 10,
                            "orthorhombic body-centered": 11,
                            "monoclinic p": 12,
                            "monoclinic base-centered": 13,
                            "triclinic p": 14}


def str2bool(v):
    """
    convert a string that has logical meaning into a boolean
    """
    return v.lower() in ("true", "t", "yes", "1", "ok")

def xml_openroot(xmlFile):
    """
    open the xml file and return the root node
    """
    root = ET.parse(xmlFile).getroot()
    return root
    
    
def xml_index(element):
    """
    print all the child found in element with their attributes (not recursively)
    """     
    for child in element:
        print child.tag, child.attrib


def get_qexml_version(tree, mydict):
    """
    get the xml format 
    """
    element = tree.find("HEADER")
    mydict['qexml_version'] = element.find("FORMAT").attrib["VERSION"].text.strip()
    mydict['qexml_after_1.4.0'] = mydict['qexml_version'] >= "1.4.0"
    

def get_eigenvalues_for_kpoint_xml(xmlfile):
    """
    return the eigenvalues associated with kpoints
    """
    kproot = ET.parse(xmlfile).getroot()
    nbnd = int(kproot.find('INFO').attrib['nbnd'])

    eigenvalues = [float(a) for a in kproot.find('EIGENVALUES').text.split('\n') if not (a.isspace() or not a)]
    occupations =[float(a) for a in kproot.find('OCCUPATIONS').text.split('\n') if not (a.isspace() or not a)]

    return eigenvalues, occupations
    

def get_k_points_xml(tree, mydict):
    """
    return the kpoints information, eigenvalues and weights
    """

    if mydict.has_key('got_kpoints'):
        if mydict['got_kpoints']:   
            return
    

    mydict['nkstot'] = int(tree.find("BRILLOUIN_ZONE").find("NUMBER_OF_K-POINTS").text.strip())
    nkstot = mydict['nkstot']
    if mydict['lsda']:
        mydict['nkpoints'] = mydict['nkstot'] / 2
    else:
        mydict['nkpoints'] = mydict['nkstot']
    
    element = tree.find("EIGENVALUES")

    mydict['kpoints_coordinates'] = []
    mydict['kpoints_weights'] = []
    mydict['kpoints_eigenvalues'] = []
    mydict['kpoints_occupations'] = []
    
    for kp in range(0, nkstot):
        tag = "./K-POINT."+str(kp+1)
        kpoint = element.find(tag)
        
        mydict['kpoints_coordinates'].append([float(a) for a in kpoint.find('K-POINT_COORDS').text.split()])
        mydict['kpoints_weights'].append(float(kpoint.find('WEIGHT').text.strip()))
        kpoints_datafile = kpoint.find('DATAFILE').attrib['iotk_link']
        eigenv, occup = get_eigenvalues_for_kpoint_xml(kpoints_datafile)
        mydict['kpoints_eigenvalues'].append(eigenv)
        mydict['kpoints_occupations'].append(occup)

        mydict['got_kpoints'] = True

        
def get_cell_xml(tree, mydict):
    """
    read in cell information: lattice vectos, cell factor..
    it adds the following keys to the dictionary mydict:
    
    NP_cell_correction, bravais_lattice, bravais_index, alat,
    lattive_vector_units, lattice_vector, lattice_rvector,
    moving_cell, cell_factor
    """
    if mydict.has_key('got_cell'):
        if mydict['got_cell']:
            return
    
    element = tree.find("CELL")
    mydict['NP_cell_correction'] = element.find("NON-PERIODIC_CELL_CORRECTION").text.strip()
    
    mydict['bravais_lattice'] = element.find("BRAVAIS_LATTICE").text.strip().lower()
    mydict['bravais_index'] = bravais_lattice_to_ibrav[mydict['bravais_lattice']]

    mydict['alat'] = float(element.find("LATTICE_PARAMETER").text.strip())

    sub_el = element.find("DIRECT_LATTICE_VECTORS")
    mydict['lattice_vector_units'] = sub_el.find("UNITS_FOR_DIRECT_LATTICE_VECTORS").attrib['UNITS'].lower().strip()
    v1 = [float(a) for a in sub_el.find("a1").text.strip().split()]
    v2 = [float(a) for a in sub_el.find("a2").text.strip().split()]
    v3 = [float(a) for a in sub_el.find("a3").text.strip().split()]
    mydict['lattice_vector'] = [v1, v2, v3]

    sub_el = element.find("RECIPROCAL_LATTICE_VECTORS")
    lattice_vector_units = sub_el.find("UNITS_FOR_RECIPROCAL_LATTICE_VECTORS").attrib['UNITS'].lower().strip()
    v1 = [float(a) for a in sub_el.find("b1").text.strip().split()]
    v2 = [float(a) for a in sub_el.find("b2").text.strip().split()]
    v3 = [float(a) for a in sub_el.find("b3").text.strip().split()]
    mydict['lattice_rvector'] = [v1, v2, v3]
                                                       

    subel = element.find("MOVING_CELL")
    if not subel is None:
        mydict['moving_cell'] = subel.text
        mydict['cell_factor'] = float(moving_cell.strip())
    else:
        mydict['moving_cell'] = None
        mydict['cell_factor'] = None
        

    mydict['got_cell'] = True

    
def get_planewaves_xml(tree, mydict):
    """
    read in plain waves associated quantities
    adds the following keys to the dictionary mydict:

    cutoff_units, cutpff_wfc, cutoff_rho, npwx (=max_number_og_gk_vectors),
    gamma_only, dfftp, ngm_g, dffts, ngms_g, got_pwaves
    
    """
    if mydict.has_key('got_pwaves'):
        if mydict['got_pwaves']:
            return
    
    element = tree.find("PLANE_WAVES")
    mydict['cutoff_units'] = element.find("UNITS_FOR_CUTOFF").text.strip().lower()
        
    mydict['cutoff_wfc'] = float(element.find("WFC_CUTOFF").text.strip())
    mydict['cutoff_rho'] = float(element.find("RHO_CUTOFF").text.strip())
    
    mydict['npwx'] = int(element.find("MAX_NUMBER_OF_GK-VECTORS").text.strip())

    mydict['gamma_only'] = str2bool(element.find("GAMMA_ONLY").text.strip())

    attrib = element.find("FFT_GRID").attrib
    mydict['dfftp'] = [float(a) for a in attrib.values()]

    mydict['ngm_g'] = int(element.find("GVECT_NUMBER").text.strip())

    attrib = element.find("SMOOTH_FFT_GRID").attrib
    mydict['dffts'] = [float(a) for a in attrib.values()]

    mydict['ngms_g'] = int(element.find("SMOOTH_GVECT_NUMBER").text.strip())

    mydict['got_pwaves'] = True

    
def get_dimensions_xml(tree, mydict):
    """
    read in "dimensions", that is variables that defines the arrays dimension
    add the following keys to the dictionary

    nat (int): num. of atoms
    nspecies (int): num. of species
    nsym (int): num. of symmetries
    nsym_brvais (int)
    lsda (bool)
    noncolin (bool)
    smearing (string)
    ntetra (int) : number of tetrahedra
    nkpstot Int) : number of kpoints (total: if lsda is true it is multiplied by 2)
    nelec (int) : num. of electrons
    nbnd (nt) : number of bands
    
    """
    if mydict.has_key('got_dimensions'):
        if mydict['got_dimensions']:
            return

    
    # IONS
    element = tree.find("IONS")
    mydict['nat'] = int(element.find("NUMBER_OF_ATOMS").text.strip())
    mydict['nspecies'] = int(element.find("NUMBER_OF_SPECIES").text.strip())
        
    # SYMMETRIES
    element = tree.find("SYMMETRIES")
        
    mydict['nsym'] = int(element.find("NUMBER_OF_SYMMETRIES").text.strip())
    mydict['nsym_bravais'] = int(element.find("NUMBER_OF_BRAVAIS_SYMMETRIES").text.strip())
    
    # PLANE WAVES
        
    # SPIN
    element = tree.find("SPIN")
    mydict['lsda'] = str2bool(element.find("LSDA").text.strip())
    mydict['noncolin'] = str2bool(element.find("NON-COLINEAR_CALCULATION").text.strip())

    # OCCUPATIONS
    element = tree.find("OCCUPATIONS")
    mydict['smearing'] = str2bool(element.find("SMEARING_METHOD").text.strip())
    el = element.find("NUMBER_OF_TETRAHEDRA")
    if el is None:
        mydict['ntetra'] = 0
    else:
        mydict['ntetra'] = int(el.text.strip())

    # BRILLOUIN ZONE
    element = tree.find("BRILLOUIN_ZONE")
    mydict['nkpstot'] = int(element.find("NUMBER_OF_K-POINTS").text.strip())
    if mydict['lsda'] is True:
        mydict['nkpstot'] *= 2
            
    # BAND STRUCTURE
    element = tree.find("BAND_STRUCTURE_INFO")
    mydict['nelec'] = float(element.find("NUMBER_OF_ELECTRONS").text.strip())
    mydict['nbnd'] = int(element.find("NUMBER_OF_BANDS").text.strip())

    mydict['got_dimensions'] = True

    
def get_ions_xml(tree, mydict):
    """
    read in quantitiesd related to the atoms
    adds the following keys to the dictionary:

    nat (int) : num of atoms  (also added by get_dimensions_xml
    nsp (int) : num opf species (the same than nspecies)

    atoms, atoms_index, atoms_tau,
    atoms_if_pos : these are lists
                   'atoms' entries are dictionaries with keys:
                   {`specie','mass', 'PPfile'}
                   'atoms_tau' and atoms_if_pos entries are lists
                   with 3 elements
    
    """

    if mydict.has_key('got_ions'):
        if mydict['got_ions']:
            return

    if not mydict['got_cell']:
        get_cell_xml(tree, mydict)

    element = tree.find("IONS")

    mydict['nat'] = int(element.find("NUMBER_OF_ATOMS").text.strip())
    mydict['nsp'] = int(element.find("NUMBER_OF_SPECIES").text.strip())

    mydict['atoms'] = []
    mydict['atoms_index'] = []
    mydict['atoms_tau'] = []
    mydict['atoms_if_pos'] = []
    

    if mydict['qexml_after_1.4.0']:
        for i in range(nsp):
            this_atom = {}
            subel = element.find("SPECIE."+str(i+1))
            
            this_atom['specie'] = subel.find("ATOM_TYPE").text.strip()
            this_atom['mass'] = float(subel.find("MASS").text.strip())
            this_atom['PPfile'] = subel.find("PSEUDO").text.strip()
            mydict['atoms'].append(this_atom)
    else:
        ValueError, "qexml version is before 1.4.0. This has not been implemented"
        
    mydict['pseudo_dir'] = element.find("PSEUDO_DIR").tex.strip()

    alat = mydict['alat']
    for i in range(nat):
        subel = element.find("ATOM."+str(i+1))

        mydict['atoms_index'].append(int(subel.attrib("INDEX")))
        mydict['atoms_tau'].append([float(s)/alat for s in subel.attrib("tau").split()])
        mydict['atoms_if_pos'].append([int(s) for s in subel.attrib("tau").split()])

    mydict['got_ions'] = True

def get_efield_xml(tree, mydict):
    """
    read electric field data
    """

    if mydict.has_key('got_efield'):
        if mydict['got_efield']:
            return

    element = tree.find("ELECTRIC_FIELD")
    if not element is None:

        mydict['has_efield'] = str2bool(element.find("HAS_ELECTRIC_FIELD").text.strip())
        mydict['has_dipole_correction'] = str2bool(element.find("HAS_DIPOLE_CORRECTION").text.strip())
        mydict['efield_dir'] = int(element.find("FIELD_DIRECTION").text.strip())
        mydict['efield_maxpos'] = float(element.find("MAXIMUM_POSITION").text.strip())
        mydict['efield_invregion'] = float(element.find("INVERSE_REGION").text.strip())
        mydict['efield_amplitude'] = float(element.find("FIELD_AMPLITUDE").text.strip())
        
    else:

        mydict['has_efield'] = False
        mydict['has_dipole_correction'] = False

    mydict['got_efield'] = True

    
def get_spin_xml(tree, mydict):
    """ """

    if mydict.has_key('got_spin'):
        if mydict['got_spin']:
            return

    element = tree.find("SPIN")
    mydict['lsda'] = str2bool(element.find("LSDA").text.strip())
    subel = element.find("NON-COLINEAR_CALCULATION")
    if not subel is None:
        mydict['noncolin'] = str2bool(subel.text.strip())
    else:
        mydict['noncolin'] = False

    if mydict['lsda']:
        spin = 2
    elif mydict['noncolin']:
        spin = 4
    else:
        spin = 1
    mydict['spin'] = spin

    if mydict['noncolin']:
        mydict['npol'] = int(element.find("SPINOR_DIM").text.strip())
    else:
        mydict['npol'] = 1
    
    subel = element.find("SPIN-ORBIT_CALCULATION")
    if not subel is None:
        mydict['spin-orbit_calculation'] = str2bool(subel.text.strip())
    else:
        mydict['spin-orbit_calculation'] = False

    subel = element.find("SPIN-ORBIT_DOMAG")
    if not subel is None:
        mydict['spin-orbit_domag'] = str2bool(subel.text.strip())
    else:
        mydict['spin-orbit_domag'] = False
        
    mydict['got_spin'] = True
    

def get_magnetization_xml(tree, mydict):
    """ """

    if mydict.has_key('got_magnetization'):
        if mydict['got_magnetization']:
            return

    element = tree.find("MAGNETIZATION_INIT")

    mydict['constraint_mag'] = int(element.find("CONSTRAINT_MAG").text.strip())
    i_cons = mydict['constraint_mag']
    nsp = int(element.find("NUMBER_OF_SPECIES").text.strip())

    mydict['mag_data'] = []
    mydict['mag_cons'] = []
    
    for i in range(nsp):
        specie = element.find("SPECIE."+str(i+1))

        mag = {}
        mag['starting_mag'] = float(specie.find("STARTING_MAGNETIZATION").text.strip())
        mag['angle1'] = float(specie.find("ANGLE1").text.strip()) * (3.141516/180.0)
        mag['angle2'] = float(specie.find("ANGLE2").text.strip()) * (3.141516/180.0)

        if (i_cons == 1) or (i_cons == 2):
            v = []
            v.append(float(specie.find("CONSTRAINT_1").tex.strip()))
            v.append(float(specie.find("CONSTRAINT_2").tex.strip()))
            v.append(float(specie.find("CONSTRAINT_3").tex.strip()))
        else:
            v = None

        mydict['mag_data'].append(mag)
        if not v is None:
            mydict['mag_cons'].append(v)

    if i_cons == 3:
        v = []
        v.append(float(specie.find("FIXED_MAGNETIZATION_1").tex.strip()))
        v.append(float(specie.find("FIXED_MAGNETIZATION_2").tex.strip()))
        v.append(float(specie.find("FIXED_MAGNETIZATION_3").tex.strip()))
        mydict['mag_cons'].append(v)

    if i_cons == 4:
        v = []
        v.append(float(specie.find("MAGNETIC_FIELD_1").tex.strip()))
        v.append(float(specie.find("MAGNETIC_FIELD_2").tex.strip()))
        v.append(float(specie.find("MAGNETIC_FIELD_3").tex.strip()))
        mydict['bfield'] = v
    else:
        mydict['bfield'] = None

    mydict['two_fermi_energies'] = str2bool(element.find("TWO_FERMI_ENERGIES").text.strip())

    if mydict['two_fermi_energies']:
        subel = element.find("TWO_FERMI_ENERGIES")
        v = float(subel.find("FIXED_MAGNETIZATION").text.strip())
        mydict['mag_cons'][0][2] = v
        mydict['nelup'] = float(subel.find("ELECTRONS_UP").text.strip())
        mydict['neldw'] = float(subel.find("ELECTRONS_DOWN").text.strip())
        mydict['ef_up'] = float(subel.find("FERMI_ENERGY_UP").text.strip())
        mydict['ef_dw'] = float(subel.find("FERMI_ENERGY_DOWN").text.strip())
        
        mydict['ef_up'] = mydict['ef_up'] * 2.0
        mydict['ef_dw'] = mydict['ef_dw'] * 2.0

    if i_cons > 0:
        mydict['lambda'] = float(element.find("LAMBDA").text.strip())
    
    mydict['got_magnetization'] = True


def get_xc_xml(tree, mydict):
    """ Exchange-Correlation """

    if mydict.has_key('got_xc'):
        if mydict['got_xc']:
            return

    if not mydict['got_ions']:
        get_ions_xml(tree, mydict)

    element = tree.find("EXCHANGE_CORRELATION")

    mydict['DFT'] = element.find("DFT").text
    subel = element.find("LDA_PLUS_U_CALCULATION")
    
    if not subel is None:
        
        nsp = int(element.find("NUMBER_OF_SPECIES").text.strip())
        mydict['lda_plus_u_kind'] = float(element.find("LDA_PLUS_U_KIND").text.strip())
        mydict['hubbard_lmax'] = float(element.find("HUBBARD_LMAX").text.strip())
        mydict['hubbard_l'] = [float(x) for x in element.find("HUBBARD_L").text.split()]
        mydict['hubbard_u'] = [float(x) for x in element.find("HUBBARD_U").text.split()]
        mydict['hubbard_j'] = [float(x) for x in element.find("HUBBARD_J").text.split()]
        mydict['hubbard_alpha'] = [float(x) for x in element.find("HUBBARD_ALPHA").text.split()]

    subel = element.find("NON_LOCAL_DF")
    if not subel is None:
        mydict['non-local_df'] = int(subel.text.strip())
        if mydict['non-local_df']==1 or mydict['non-local_df']== 2:
            mydict['non-local_df'] = element.find("VDW_KERNEL_NAME").text
    else:
        mydict['non-local_df'] = 0
        
    mydict['got_xc'] = True
    

def get_brillouin_zone_xml(tree, mydict):
    """ """

    if mydict.has_key('got_broullouin'):
        if mydict['got_brillouin']:
            return

    element = tree.find("BRILLOUIN_ZONE")
    num_kpts = int(element.find("NUMBER_OF_K-POINTS").text.strip())
    nkstot = num_kpts

    if mydict['lsda']:
        nkstot = nkstot * 2

    mydict['nk'] = [int(x) for x in element.find("MONKHORST_PACK_GRID").attr.values()]
    mydict['k'] = [int(x) for x in element.find("MONKHORST_PACK_OFFSET").attr.values()]

    mydict['xk'] = []
    mydict['wk'] = []

    v = []
    w = []
    for i in range(num_kpts):
        element.find("K-POINT."+str(i+1)).attr
        v.append([float(x) for x in attr['XYZ'].split()])
        w.append(float(attr['WEIGHT']))

    mydict['xk'] = v.copy()
    mydict['wk'] = w.copy()
    
    if mydict['lsda']:
        mydict['xk'].append(v.copy())
        mydict['wk'].append(w.copy())


    num_kpts_start = int(element.find("STARTING_K-POINTS").text.strip())        
    mydict['xk_start'] = []
    mydict['wk_start'] = []
    for i in range(num_kpts_start):
        element.find("K-POINT_START."+str(i+1)).attr
        mydict['xk_start'].append([float(x) for x in attr['XYZ'].split()])
        mydict['wk_start'].append(float(attr['WEIGHT']))

    # for completeness, Bravais-symmetries should also be read here // look in pw_restart.f90 from PW source code dir
        
    mydict['got_brillouin'] = True


def get_occupations_xml(tree, mydict):
    """ """

    if mydict.has_key('got_occupations'):
        if mydict['got_occupations']:
            return

    element = tree.find("OCCUPATIONS")

    subel = element.find("SMEARING_METHOD")
    if not subel is None:
        mydict['lgauss'] = str2bool(subel.text.strip())
    else:
        mydict['lgauss'] = False
    
    if mydict['lgauss']:
        mydict['ngauss'] = int(element.find("SMEARING_TYPE").text.strip())
        if mydict['ngauss'] == 0:
            mydict['smearing'] = "gaussian"
        if mydict['ngauss'] == 1:
            mydict['smearing'] = "Methfessel-Paxton"
        if mydict['ngauss'] == 2:
            mydict['smearing'] = "Marzari-Vanderbilt"
        if mydict['ngauss'] == 3:
            mydict['smearing'] = "Fermi-Dirac"
            
        mydict['degauss'] = float(element.find("SMEARING_PARAMETER").text.strip())
        mydict['degauss'] = mydict['degauss'] * 2.0
        
    else:
        mydict['ngauss'] = 0
        mydict['degauss'] = 0.0

    mydict['ltetra'] = str2bool(element.find("TETRAHEDRON_METHOD").text.strip())
    if mydict['ltetra']:
        mydict['ntetra'] = int(element.find("NUMBER_OF_TETRAHEDRA").text.strip())
        ntetra = mydict['ntetra']
        mydict['tetrahedra'] = []
        for i in range (ntetra):
             mydict['tetrahedra'].append([float(x) for x in element.find("TETRAHEDRON"+str(i+1)).text.split()])
    else:
        mydict['ntetra'] = 0

    subel = element.find("FIXED_OCCUPATIONS")
    if not subel is None:
        mydict['tfixed_occ'] = str2bool(subel.text.strip())
    else:
        mydict['tfixed_occ'] = False

    if mydict['tfixed_occ']:
        subel = element.find("INFO")
        if not subel is None:
            mydict['nupdwn'] = [subel.attr['nstates_up'], subel.attr['nstates_down']]
            
        else:
            mydict['nupdwn'] = [mydict['nbnd'], mydict['nbnd']]

        mydict['finp'] = [float(x) for x in element.find("INPUT_OCC_UP").text.split()]
        if mydict['lsda']:
            mydict['finp'].append([float(x) for x in element.find("INPUT_OCC_DOWN").text.split()])

        
    mydict['ltetra'] = str2bool(element.find("TETRAHEDRON_METHOD").text.strip())
        
    mydict['got_occupations'] = True


def get_band_structure_xml(tree, mydict):
    """ """

    if mydict.has_key('got_band_structure'):
        if mydict['got_band_structure']:
            return

    if not mydict['got_spin']:
        get_spin_xml(tree, mydict)
    if not mydict['got_brillouin']:
        get_brillouin_zone_xml(tree, mydict)
    
    info = tree.find('BAND_STRUCTURE_INFO')
    
    mydict['nelec'] = int(element.find("NUMBER_OF_ELECTRONS").text.strip())
    subel = element.find("NUMBER_OF_ATOMIC_WFC")
    if not subel is None:
        mydict['natomwfc'] = int(element.find("NUMBER_OF_ATOMIC_WFC").text.strip())
    else:
        mydict['natomwfc'] = 0

    mydict['nbnd'] = int(element.find("NUMBER_OF_BANDS").text.strip())

    subel = element.find("FERMI_ENERGY")
    if not subel is None:
        mydict['ef'] = float(subel.text.strip()) * 2
    else:
        mydict['ef'] = 0

    subel = element.find("TWO_FERMI_ENERGIES")
    if not subel is None:
        two_ef = str2bool(subel.text.strip())
    else:
        two_ef = False

    if two_ef:
        mydict['ef_up'] = float(elememt.find("FERMI_ENERGY_UP").text.strip()) * 2
        mydict['ef_dw'] = float(elememt.find("FERMI_ENERGY_DOWN").text.strip()) * 2
        
    get_k_points_xml(tree, mydict)

    mydict['got_band_structure'] = True
    mydict['got_ef'] = True

    
def get_fermi_energy_xml(tree, mydict):
    """
    read in the fermi energy
    """

    if mydict.has_key('got_ef'):
        if mydict['got_ef']:
            return

    element = tree.find("BAND_STRUCTURE_INFO")
    if not element.find("FERMI_ENERGY") is None:
        mydict['ef'] = float(element.find("FERMI_ENERGY").text.strip()) * 2
    else:
        mydict['ef'] = 0.0
        
    mydict['got_ef'] = True


# ......................................................................................
    
def read_quantumespresso_xmloutput(filename, action, mydict=None):
    """
    read the output text file from a QE run and extracts physical quantities
    arguments:
      - filename : that is the name of the xml output file you want to read in
      - action : string: that specifies what exactly you want to read in. It basically
                         follows the convention from pw_restart.f90
                         'all' : read in all the data
                         'ef'  : read in just the fermi energy
                         'kpts': read in only the k-points with eigenvalues and weights
                         'config': read in the cell and atom parameters
                         ... : just check the lines below for further information
                         
    the argument 
    returns a dictionary that contains the results.
    """

    calls = {}
    call['dim'] = get_dimensions_xml
    call['cell'] = get_cell_xml
    call['pw'] = get_planewaves_xml
    call['ions'] = get_ions_xml
    call['spin'] = get_spin_xml
    call['init_mag'] = get_magnetization_xml
    call['xc'] = get_xc_xml
    call['occ'] = None
    calls['kpts'] = get_k_points_xml
    call['bz'] = get_brillouin_zone_xml
    call['bs'] = get_band_structure_xml
    call['wfc'] = get_wfc_xml
    call['efield'] = get_efield_xml
    call['symm'] = None
    call['rho'] = None
    call['pseudo'] = get_ions_xml
    call['occupations'] = get_occupations_xml

    actions = {}
    actions['dim'] = ['dim', 'bz']
    actions['pseudo'] = ['ions']
    actions['config'] = ['cell, ions']
    actions['rho'] = None
    actions['wave'] = ['pw', 'wfc']
    actions['nowave'] = ['cell', 'pw','ions', 'spin', 'init_mag', 'xc', 'occ', 'bz', 'bs', 'symm', 'efield', 'occupations']
    actions['all'] = actions['nowave'].append('rho')
    actions['ef'] = ['ef']
    actions['kpts'] = ['kpts']

    if mydict is None:
        mydict = {}
    elif action == 'all':
        mydict.clear()
    
    if action in actions.keys():
        
        root = xml_openroot(filename)
        get_dimensions_xml(root, mydict)
        
        for A in actions[action]:
            calls[A](root, mydict)

    else:
        ValueError, "the action you specified ("+action+") is not present in my archive"

    return mydict



# -------------------------------------------------------------------------------------------------------------------------------
#
#  read infile
#
# -------------------------------------------------------------------------------------------------------------------------------    


# ......................................................................................

def write_atoms(infilename, myprefix="dry_pwscf", verbose = False):
    """
    this function builds an Atoms object in ASE starting
    from a .in file, using the data collected by
    read_in_file()
    blocks must be the dictionary returned by read_in_file()

    it returns an Atoms object
    """

    new_infile = infilename+'.dry'
    
    blocks = read_in_file(infilename)
    # changes the prefix in myprefix
    if blocks.has_key('prefix'):
        _command_string = "sed 's|prefix[ ]*=[ ]*[a-z/0-9]*|prefix = "+myprefix+",|' < "+infilename+" > "+new_infile
    else:
        _command_string = "sed 's|\&control|\&control\\n    prefix = "+myprefix+",|' < "+infilename+" > "+new_infile

    print _command_string
    os.system(_command_string)
        
    
    if blocks['system']['ibrav'] == 0:
        if verbose:
            print "\t[write_atoms] ibrav = 0 found\n"
        
        if blocks['system'].has_key('alat'):
            alat = blocks['system']['alat']
        elif blocks['system'].has_key('celldm(1)'):
            alat = blocks['system']['celldm(1)']
        elif blocks['system'].has_key('a'):
            alat = blocks['system']['a']
        else:
            raise ValueError, "something got wrong: it seems that ibrav=0 but neither 'alat' nor 'celldm(1)' nor 'a' are present"

        # creates the cell as an array
        cell = np.array(blocks['CELL_PARAMETERS']['cell'])
        if verbose:
            print "\t[write_atoms] found cell:\n"
            for R in cell:
                print "\t\t %f %f %f" % tuple(R)
            
        # convert the cell in bohr units
        # (note that alat or celldm(1) are supposed to be already in bohr)
        if blocks['CELL_PARAMETERS']['attrib'] == 'alat' or blocks['CELL_PARAMETERS']['attrib'] == None:
            cell = cell * alat
        elif blocks['CELL_PARAMETERS']['attrib'] == 'angstrom':
            cell = cell * convert('a','bohr')
        if verbose:
            print "\t[write_atoms] cell rescaled is:\n"
            for R in cell:
                print "\t\t %f %f %f" % tuple(R)
            
            
    else:
        # it is needed to start qe in dry mode to obtains the
        # cell parameters, because ibrav >0
        # QE will do this for us

        if verbose:
            print "\t[write_atoms] ibrav > 0 found. Invoking QE\n"
        
        if 'ESPRESSO_ROOT' in os.environ:
                rootdir = os.environ['ESPRESSO_ROOT']
                bindir = rootdir + '/PW/src/'
        else:
            rootdir= '$HOME/espresso/'
            bindir = rootdir + '/PW/src/'

        execbin = "pw.x"
        if not os.access(bindir+execbin, os.F_OK):
            raise IOError, "binary %s does not exist" % (bindir+execbin)
        if not os.access(bindir+execbin, os.X_OK):
            raise IOError, "you do not have execution permission on the binary %s" % bindir+execbin

        # run a dry run: only the xml data file is written, with the information about the cell
        tempfilename = myprefix+".EXIT"
        fh_temp = open(tempfilename, "w")
        fh_temp.close()
        _commandstring = '%s < %s > %s' % ( bindir+execbin, new_infile, new_infile+'.out')
        exitcode = os.system(_commandstring)
        
        # read in cell information
        if verbose:
            print "\t[write_atoms] now read-in the xml output\n"
        root = ET.ElementTree(file=myprefix+'.save/data-file.xml').getroot()
        mydict = {}
        get_cell_xml(root, mydict)

        alat = mydict['alat']
        # the cell is already in bohr
        cell = mydict['lattice_vector']

        if verbose:
            print "\t[write_atoms] found cell:\n"
            for R in cell:
                print "\t\t %f %f %f" % tuple(R)

    # now actually creates the Atoms structure
        
    atoms = Atoms()
    for A in blocks['ATOMIC_POSITIONS']['atoms']:
        print "\t[write_atoms] adding atom: %s\n" % A['symbol']
        atoms.append(Atom(A['symbol'], tuple(A['pos'])))
    print "\t[write_atoms] attaching cell:\n"
    for R in cell:
        print "\t\t %f %f %f" % tuple(R)
    
    atoms.set_cell(cell, scale_atoms=True)

    return atoms


# ......................................................................................

def read_in_file(infilename):
    """
    this function read a quantum espresso  .in file
    and returns a dictionary for all the namelists
    and the cards found.
    it correctly deals lines with more than one
    comma-separated parameters

    it returns a dictionary, whos keys are the name
    of the namelists and cards.
    In turn, each value is a dictionary with all the
    values found.

    as for the namecards:
    (+) atomic_species is a dictionary whos keys are progeressive
        integers and whose values are dictionaries with the following
        keys:
        'symbol', 'attrib', 'mass', 'PP'
        example:
        {'1': {'symbol: 'H', 'attrib': 1, 'mass': 1.0008, 'PP': "pseudo_file_for_H"},
        '2': {'symbol: 'C', 'attrib': None, 'mass': 1, 'PP': "pseudo_file_for_C"},
        '3': {'symbol: 'H', 'attrib': 2, 'mass': 1.0008, 'PP': "different_pseudo_file_for_H"},
        ...}
        the key 'attrib' is set in case there is more than 1 entry for a specie, for instance
        to specify different PP:
        ATOMIC_SPECIES
         H1 1.0008 pseudo_file_for_H
         H2 1.0008 different_pseudo_file_for_H
         C  12 pseudo_file_for_C

    (+) atomic_positions is a dictionary whos keys are progressive
    integers and whose values are dictionaries with the following keys:

    example: 'symbol', 'pos', 'if_pos'
    {'1': {'symbol: 'H', 'pos': [1.0, 1.0, 1.0], 'if_pos': [0.1, 0.1, 0.1]},
     ...}
     the key 'if_pos' may have value None if it wasn't present

     (+) all the other namecards are dictionaries with a unique key, 'elements',
        that is a list of all the lines inside that card
     
    """

    if not os.access(infilename, os.F_OK):
        raise IOError, "the .in file %s is not present\n" % infilename

    fh = open(infilename)
    alltext = fh.read()
    alllines = alltext.split('\n')
    
    card_labels = ['ATOMIC_SPECIES',
                   'ATOMIC_POSITIONS',
                   'CELL_PARAMETERS',
                   'K_POINTS',
                   'CONSTRAINTS'
                   'OCCUPATIONS']

    blocks = {}
    blocks['sparse']=[]
    isnamelist = False
    iscard = False

    for line in alllines:

        # determines whether we are inside a namelist or a card
        if line.strip().startswith('&'):
            isnamelist = True
            blockname = line.strip()[1:]
            blocks[blockname] = {}

        elif line.strip().startswith('/'):
            isnamelist = False

        elif line.strip():
            if line.split()[0].strip() in card_labels:
                iscard = True
                
                blockname = line.strip().split()[0]
                blocks[blockname] = {}

                try:
                    attrib = line[line.find(line.strip().split()[1]):].strip()
                    if attrib.startswith('{'):
                        attrib = attrib[1:]
                    if attrib.endswith('}'):
                        attrib = attrib[:-1]
                except:
                    attrib = None
                blocks[blockname]["attrib"] = attrib

            else:     # if in a namelist, isolate keywords and values

                if isnamelist:
                    tokens = line.split(',')
                    for t in tokens:
                        if t.strip():

                            key, svalue = t.strip().split('=')

                            key = key.strip()
                            svalue = svalue.strip()
                            if svalue.endswith(','):
                                value = svalue[:-1].strip()
                            try:
                                value=int(svalue)
                            except ValueError:
                                try:
                                    value=float(svalue)
                                except ValueError:
                                    if svalue.lower() in ['.true.', 'true' , 't']:
                                        value = True
                                    elif svalue.lower() in ['.false.', 'false' , 'f']:
                                        value = False
                                    else:
                                        value = str(svalue)

                            blocks[blockname][key] = value

                elif iscard:

                    if blockname == "ATOMIC_SPECIES":                        # -- SPECIES
                        tokens = line.split()
                        if(blocks[blockname].has_key('count')):
                            blocks[blockname]['count'] = blocks[blockname]['count']+1
                        else:
                            blocks[blockname]['count'] = 1
                            
                        symbol = tokens[0].strip()
                        
                        if symbol.find('-')> 0 or symbol.find('_')>0:
                            ssymbol = symbol.split('-')
                            if len(ssymbol) == 1:
                                ssymbol = symbol.split('_')
                            symbol=ssymbol[0]
                            symbol_attrib=ssymbol[1]
                        elif not symbol.isalpha():
                            digits="0123456789"
                            found=min([symbol.index(d) for d in digits if d in symbol])
                            symbol_attrib = int(symbol[found:])
                            symbol = symbol[:found]
                        else:
                            symbol_attrib = None
                            
                        blocks[blockname][str(blocks[blockname]['count'])] = {'symbol': symbol,
                                                                              'attrib': symbol_attrib,
                                                                               'mass' : float(tokens[1]),
                                                                               'PP' : tokens[2].strip()}

                    elif blockname == "ATOMIC_POSITIONS":                    # -- POSITIONS
                        tokens = line.split()
                        if(blocks[blockname].has_key('count')):
                            blocks[blockname]['count'] = blocks[blockname]['count']+1
                        else:
                            blocks[blockname]['count'] = 1
                            blocks[blockname]['atoms'] = []

                        idx = blocks[blockname]['count']
                        if len(tokens) > 4:
                            if_pos = [float(eval(s)) for s in tokens[4:]]
                        else:
                            if_pos = None

                        blocks[blockname]['atoms'].append({'symbol': tokens[0].strip(),
                                                          'pos': [float(eval(s)) for s in tokens[1:4]],
                                                          'if_pos': if_pos})

                    elif blockname == "CELL_PARAMETERS":                     # -- CELL PARAMETERS

                        tokens = line.split()

                        if not blocks[blockname].has_key('cell'):
                            blocks[blockname]['cell'] = []

                        blocks[blockname]['cell'].append([float(s) for s in tokens[0:3]])
                        if(len(blocks[blockname]['cell'])) == 3:
                            iscard = False                                                

                    else:
                        if not blocks[blockname].has_key('elements'):
                            blocks[blockname]['elements'] = []
                        blocks[blockname]['elements'].append(line)

                else:
                    blocks['sparse'].append(line)


    fh.close()
    return blocks


def read_log(log_filename):
    """Read a log file as produced from a QuantenEspresso calculation
    and return an Atoms object.

    """
    calc = calculator()
    p = os.popen('grep -n Giannozzi {log_filename} 2>/dev/null | tail -1'.format(**locals()), 'r')

    try:
        lline = p.readline()
        print(r'{log_filename} preadline {lline}'.format(**locals()))
        n = int(lline.split()[0].strip(':'))
    except:
        print >>sys.stderr, 'No valid pw-log at {log_filename} found.'.format(**locals())
        p.close()
        return
    p.close()

    s = open(log_filename, 'r')
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
    i = 0
    while a[:12] != '     site n.':
        a = s.readline()
        if not a:
            print >>sys.stderr, 'no regular beginning of calculation found'
            break
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
    atoms0 = Atoms(syms, pos, cell=cell * alat, pbc=(1, 1, 1))

    atoms0.get_calculator = calc.get_calculator
    atoms0.get_potential_energy = calc.get_potential_energy
    atoms0.get_forces = calc.notimpl
    atoms0.get_stress = calc.notimpl
    atoms0.get_charges = calc.notimpl

    # get total energy at first ionic step
    en = get_total_energy(s)
    if en is not None:
        calc.set_energy(en)
    else:
        print >>sys.stderr, 'no total energy found'
        return

    a = s.readline()

    if a != '':
        traj = [atoms0]
    else:
        return atoms0

    while a != '':
        atoms = copy.deepcopy(atoms0)
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
            # get total energy at 2nd,3rd,...nth ionic step
            en = get_total_energy(s)
            calc.set_energy(en)

            if en is not None:
                traj.append(atoms)
            else:
                break
        else:
            for i in range(3):
                cell[i][:] = s.readline().split()
            atoms.set_cell(cell * alat, scale_atoms=False)
        a = s.readline()



    return traj
