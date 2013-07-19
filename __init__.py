svnver = 'SVNVERSION'
import os

try:
    import espsite
except:
    print '*** ASE Espresso requires a site-specific espsite.py in PYTHONPATH.'
    print '*** You may use espsite.py.example in the svn checkout as a template.'
    raise ImportError
site = espsite.config()

from ase.calculators.general import Calculator
import atexit
import sys, string
import numpy as np
from types import FileType, StringType
from constants import *
from utils import *
from subdirs import *

# ase controlled pw.x's register themselves here, so they can be
# stopped automatically
espresso_calculators = []


class espresso(Calculator):
    def __init__(self,
                 pw = 350.0,
                 dw = None,
                 nbands = -10, 
                 kpts = (1,1,1),
                 kptshift = (0,0,0),
                 mode = 'ase3',
                 opt_algorithm = 'ase3',
                 fmax = 0.05,
                 cell_dynamics = None,
                 press = None, # target pressure
                 dpress = None, # convergence limit towards target pressure
                 cell_factor = None,
                 dontcalcforces = False,
                 nosym = False,
                 noinv = False,
                 nosym_evc = False,
                 no_t_rev = False,
                 xc = 'PBE',
                 psppath = None,
                 spinpol = False,
                 noncollinear = False,
                 spinorbit = False,
                 outdir = None,
                 txt = None,
                 calcstress = False,
                 smearing = 'fd',
                 sigma = 0.1,
                 fix_magmom = False,
                 U = None,
                 J = None,
                 U_alpha = None,
                 U_projection_type = 'atomic',
                 tot_charge = None, # +1 means 1 e missing, -1 means 1 extra e
                 tot_magnetization = -1, #-1 means unspecified, 'hund' means Hund's rule for each atom
                 occupations = 'smearing', # 'smearing', 'fixed', 'tetrahedra'
                 dipole = {'status':False},
                 field = {'status':False},
                 output = {'disk_io':'default',  # how often espresso writes wavefunctions to disk
                           'avoidio':False,  # will overwrite disk_io parameter if True
                           'removewf':True,
                           'wf_collect':False},
                 convergence = {'energy':1e-6,
                                'mixing':0.7,
                                'maxsteps':100,
                                'diag':'david'},
                 startingpot = None,
                 startingwfc = None,
                 parflags = None,
                 onlycreatepwinp = None, #specify filename to only create pw input
                 single_calculator = True, #if True, only one espresso job will be running
                 procrange = None, #let this espresso calculator run only on a subset of the requested cpus
                 numcalcs = None,  #used / set by multiespresso class
		 verbose = 'low'):
        
        self.outdir= outdir
        self.onlycreatepwinp = onlycreatepwinp 
        self.pw = pw
        self.dw = dw
        self.nbands = nbands
        self.kpts = kpts
        self.kptshift = kptshift
        self.calcmode = mode
        self.opt_algorithm = opt_algorithm
        self.fmax = fmax
        self.cell_dynamics = cell_dynamics
        self.press = press
        self.dpress = dpress
        self.cell_factor = cell_factor
        self.dontcalcforces = dontcalcforces
        self.nosym = nosym
        self.noinv = noinv
        self.nosym_evc = nosym_evc
        self.no_t_rev = no_t_rev
        self.xc = xc
        self.smearing = smearing
        self.sigma = sigma
        self.spinpol = spinpol
        self.noncollinear = noncollinear
        self.spinorbit = spinorbit
        self.fix_magmom = fix_magmom
        self.tot_charge = tot_charge
        self.tot_magnetization = tot_magnetization
        self.occupations = occupations
        self.outdir = outdir
        self.calcstress = calcstress
        self.psppath = psppath
        self.dipole = dipole
        self.field = field
        self.output = output
        self.convergence = convergence
        self.startingpot = startingpot
        self.startingwfc = startingwfc
        self.verbose = verbose
        self.U = U
        self.J = J
        self.U_alpha = U_alpha
        self.U_projection_type = U_projection_type
        if parflags is None:
            self.parflags = ''
        else:
            self.parflags = parflags
	self.single_calculator = single_calculator
        self.txt = txt

        self.mypath = os.path.abspath(os.path.dirname(__file__))
        self.writeversion = True

        self.atoms = None
        self.sigma_small = 1e-13
        self.started = False
        self.got_energy = False
        self.only_init = False

        # Variables that cannot be set by inputs
        self.nvalence=None
        self.nel = None 
        # Auto create variables from input
        self.input_update() 

        # Initialize lists of cpu subsets if needed
        if procrange is None:
            self.proclist = False
        else:
            self.proclist = True
            procs = site.procs + []
            procs.sort()
            nprocs = len(procs)
            self.myncpus = nprocs / numcalcs
            i1 = self.myncpus * procrange
            self.mycpus = self.localtmp+'/myprocs%04d.txt' % procrange
            f = open(self.mycpus, 'w')
            for i in range(i1,i1+self.myncpus):
                print >>f, procs[i]
            f.close()


    def input_update(self):
        """ Run initialization functions, such that this can be called if variables in espresso are
        changes using set or directly. 
        """
        self.create_outdir() # Create the tmp output folder 

        #sdir is the directory the script is run or submitted from
        self.sdir = getsubmitorcurrentdir(site)

        if self.dw is None:
            self.dw = 10. * self.pw
        else:
            assert self.dw >= self.pw

        if self.psppath is None:
            try:
                self.psppath = os.environ['ESP_PSP_PATH']
            except:
                print 'Unable to find pseudopotential path.  Consider setting ESP_PSP_PATH environment variable'
                raise
        if self.dipole is None:
            self.dipole = {'status':False}
        if self.field is None:
            self.field = {'status':False}
        
        if self.convergence is None:
            self.conv_thr = 1e-6/rydberg
        else:
            if self.convergence.has_key('energy'):
                self.conv_thr = self.convergence['energy']/rydberg
            else:
                self.conv_thr = 1e-6/rydberg
        self.started = False
        self.got_energy = False

    def create_outdir(self):
        if self.onlycreatepwinp is None:
            self.localtmp = mklocaltmp(self.outdir, site)
            if not self.txt:
                self.log = self.localtmp+'/log'
            elif self.txt[0]!='/':
                self.log = self.sdir+'/log'
            else:
                self.log = self.txt
            self.scratch = mkscratch(self.localtmp, site)
            if self.output is not None and self.output.has_key('removewf'):
                removewf = self.output['removewf']
            else:
                removewf = True
            atexit.register(cleanup, self.localtmp, self.scratch, removewf, self, site)
            self.cancalc = True
        else:
            self.pwinp = self.onlycreatepwinp
            self.localtmp=''
            self.cancalc = False


    def set(self,  **kwargs):
        """ Define settings for the Quantum Espresso calculator object after it has been initialized. 
        This is done in the following way:

        >> calc = espresso(...)
        >> atoms = set.calculator(calc)
        >> calc.set(xc='BEEF')

        NB: No input validation is made
        """
        for key, value in kwargs.items():
            if key == 'outdir':
                self.create_outdir()
                self.outdir = value
            if key == 'startingpot':
                self.startingpot = value
            if key == 'startingwfc':
                self.startingwfc = value
            if key == 'U_alpha':
                self.U_alpha = value
            if key == 'U':
                self.U = value
            if key == 'U_projection_type':
                self.U_projection_type = value
            if key == 'xc':
                self.xc = value
            if key == 'pw':
                self.pw = num2str(value)
            if key == 'dw':
                self.dw =num2str(value)
            if key == 'output':
                self.output = value
            if key == 'convergence':
                self.convergence = value
            if key == 'kpts':
                self.kpts = value
            if key == 'kshift':
                self.kshift = value
        self.input_update()
        self.recalculate = True

    def __del__(self):
        try:
            self.stop()
        except:
            pass

    def atoms2species(self):
        """ Define several properties of the quantum espresso species from the ase atoms object.
        
        Constructs a dictionary (self.specdict) with the elements and their associated properties.
        Each element has the keys: 
        'mass', 
        'magmoms' : a list with the number of unique magnetic moments, used to define # of species
        
        Constructs self.specprops which contain species labels, masses, magnetic moments, and positions.
        Also defines self.species and self.nspecies.
        """
        atypes = list( set(self.atoms.get_chemical_symbols()) )
        
        aprops = zip(self.atoms.get_chemical_symbols(), self.atoms.get_masses(),
                     self.atoms.get_initial_magnetic_moments(), 
                     self.atoms.get_scaled_positions())
        
        magmoms = self.atoms.get_initial_magnetic_moments()
        
        #### CREATE DICTIONARY FOR EACH ATOM TYPE ####
        typedict = {}
        for atype in atypes:
            typedict[atype] = {}
            maglist  = []
            masslist = []
            Ulist    = []
            Jlist    = []
            U_alphalist = []
            for i, atom in enumerate(self.atoms):
                if atom.symbol == atype:
                    maglist.append(magmoms[i])
                    masslist.append(atom.mass)
                    if self.U is not None:
                        if i in self.U:
                            Ulist.append(self.U[i])
                        elif atom.symbol in self.U:
                            Ulist.append(self.U[atom.symbol])
                        else:
                            Ulist.append(num2str(0.))
                    if self.J is not None:
                        if i in self.J:
                            Jlist.append(self.J[i])
                        elif atom.symbol in self.J:                          
                            Jlist.append(self.J[atom.symbol])
                        else:
                            Jlist.append(num2str(0.))
                    if self.U_alpha is not None:
                        if i in self.U_alpha:
                            U_alphalist.append(self.U_alpha[i])
                        elif atom.symbol in self.U_alpha:
                            U_alphalist.append(self.U_alpha[atom.symbol])
                        else:
                            U_alphalist.append(num2str(0.))
            
            typedict[atype]['magmoms'] = list( set(maglist) )
            typedict[atype]['mass'] = list( set(masslist) )
            #### uniqueness identify which atoms of same type that are unique
            un = uniqueness(maglist,masslist)
            if self.U is not None:
                typedict[atype]['U'] = list( set(Ulist))
                un = uniqueness(un,Ulist)
            if self.J is not None:
                typedict[atype]['J'] = list( set(Jlist))
                un = uniqueness(un,Jlist)
            if self.U_alpha is not None:
                typedict[atype]['U_alpha'] = list( set(U_alphalist))
                un = uniqueness(un,U_alphalist)
            typedict[atype]['indexes']=[int(kk) for kk in un]


        #### CREATE INDICES FOR EACH ATOM TYPE WITH A UNIQUE SET OF MAGMOM, U, 
        ##   J, U_alpha STARTING AT 1 ####


        speciesindex = []
        for type, info in typedict.iteritems():
            for val in info['indexes']:
                speciesindex.append(type+str(val))
        """
        for type, info in typedict.iteritems():
            for num in 
            #tcount = 1
            #for mag in info['magmoms']:
            #    speciesindex.append(type+str(tcount))
            #    tcount += 1
        """
        #### UPDATE THE SPECIES PROPERTIES TO INCLUDE THE SPECIES ID #####
        specprops = []
        index_counter = np.zeros(len(atypes)).astype(int)
        jj = 0
        for ii, a in enumerate(aprops):
            atypeindex = np.where(np.array(atypes)==a[0])[0][0]
            index = typedict[a[0]]['indexes'][index_counter[atypeindex]]
            specprops.append( (a[0]+str(index),a[1],a[2],a[3]))
            index_counter[atypeindex]+=1
        
        self.specdict  = typedict
        self.species   = np.unique(speciesindex)
        self.nspecies  = len(self.species)
        self.specprops = specprops

    def get_nvalence(self):
        nel = {}
        for x in self.species:
            el = x.strip('0123456789')
            #get number of valence electrons from pseudopotential or paw setup
            p = os.popen('egrep -i \'z\ valence|z_valence\' '+self.psppath+'/'+el+'.UPF | tr \'"\' \' \'','r')
            for y in p.readline().split():
                if y[0].isdigit() or y[0]=='.':
                    nel[el] = int(round(float(y)))
                    break
            p.close()
        nvalence = np.zeros(len(self.specprops))
        for i,x in enumerate(self.specprops):
            nvalence[i] = nel[x[0].strip('0123456789')]
        #if not self.spinpol:
        #    nvalence /= 2 
        return nvalence, nel

    def write_charge(self, f):
        if self.tot_charge != None:
           print >>f, '  tot_charge='+str(self.tot_charge)+','

    def writeinputfile(self, filename='pw.inp', mode=None,
        overridekpts=None, overridekptshift=None, suppressforcecalc=False,
        usetetrahedra=False):
        if self.atoms is None:
            raise ValueError, 'no atoms defined'
        if self.cancalc:
            fname = self.localtmp+'/'+filename
            #f = open(self.localtmp+'/pw.inp', 'w')
        else:
            fname = self.pwinp
            #f = open(self.pwinp, 'w')
        f = open(fname, 'w')
        
        ### &CONTROL ###
        if mode is None:
            if self.calcmode=='ase3':
                print >>f, '&CONTROL\n  calculation=\'relax\',\n  prefix=\'calc\','
            elif self.calcmode=='hund':
                print >>f, '&CONTROL\n  calculation=\'scf\',\n  prefix=\'calc\','
            else:
                print >>f, '&CONTROL\n  calculation=\''+self.calcmode+'\',\n  prefix=\'calc\','
        else:
                print >>f, '&CONTROL\n  calculation=\''+mode+'\',\n  prefix=\'calc\','

        print >>f, '  pseudo_dir=\''+self.psppath+'\','
        print >>f, '  outdir=\'.\','
        efield = (self.field['status']==True)
        dipfield = (self.dipole['status']==True)
        if efield or dipfield:
            print >>f, '  tefield=.true.,'
            if dipfield:
                print >>f, '  dipfield=.true.,'
        if not self.dontcalcforces and not suppressforcecalc:
            print >>f, '  tprnfor=.true.,'
            if self.calcstress:
                print >>f, '  tstress=.true.,'
            if self.output is not None:
                if self.output.has_key('avoidio'):
                    if self.output['avoidio']:
                        self.output['disk_io'] = 'none'
                if self.output.has_key('disk_io'):
                    if self.output['disk_io'] in ['high', 'low', 'none']:
                        print >>f, '  disk_io=\''+self.output['disk_io']+'\','

                if self.output.has_key('wf_collect'):
                    if self.output['wf_collect']:
                        print >>f, '  wf_collect=.true.,'
        if self.opt_algorithm!='ase3' or not self.cancalc:
            # we basically ignore convergence of total energy differences between
            # ionic steps and only consider fmax as in ase
            print >>f, '  etot_conv_thr=1d0,' 
            print >>f, '  forc_conv_thr='+num2str(self.fmax/rydberg_over_bohr)+','

        ### &SYSTEM ###
        print >>f, '/\n&SYSTEM\n  ibrav=0,\n  celldm(1)=1.8897261245650618d0,'
        print >>f, '  nat='+str(self.natoms)+','
        self.atoms2species() #self.convertmag2species()
        print >>f, '  ntyp='+str(self.nspecies)+',' #str(len(self.msym))+','
        self.write_charge(f)
        if self.calcmode!='hund':
            inimagscale = 1.0
        else:
            inimagscale = 0.9
        if self.fix_magmom:
            assert self.spinpol
            self.totmag = self.summed_magmoms
            print >>f, '  tot_magnetization='+num2str(self.totmag*inimagscale)+','
        elif self.tot_magnetization != -1:
            if self.tot_magnetization != 'hund':
                self.totmag = self.tot_magnetization
            else:
                from atomic_configs import hundmag
                self.totmag = sum([hundmag(x) for x in self.atoms.get_chemical_symbols()])
            print >>f, '  tot_magnetization='+num2str(self.totmag*inimagscale)+','
        print >>f, '  ecutwfc='+num2str(self.pw/rydberg)+','
        print >>f, '  ecutrho='+num2str(self.dw/rydberg)+','
        if self.nbands is not None:
            if self.nbands>0:
                print >>f, '  nbnd='+str(int(self.nbands))+','
            else:
                if self.nvalence == None:
                     self.nvalence, self.nel =  self.get_nvalence()
                print >>f, '  nbnd='+str(int(np.sum(self.nvalence)-self.nbands))+','
        if usetetrahedra:
            print >>f, '  occupations=\'tetrahedra\','
        else:
            if abs(self.sigma)>1e-13:
                print >>f, '  occupations=\''+self.occupations+'\','
                print >>f, '  smearing=\''+self.smearing+'\','
                print >>f, '  degauss='+num2str(self.sigma/rydberg)+','
            else:
                if self.spinpol:
                    assert self.fix_magmom
                print >>f, '  occupations=\'fixed\','
        if self.spinpol:
            print >>f, '  nspin=2,'
            spcount  = 1
            if self.nel == None:
                self.nvalence, self.nel = self.get_nvalence()
            for species in self.species: # FOLLOW SAME ORDERING ROUTINE AS FOR PSP                
                magindex = int(string.join([i for i in species if i.isdigit()],''))
                el  = species.strip('0123456789')
                mag = self.specdict[el]['magmoms'][magindex-1]/self.nel[el]
                assert np.abs(mag) <= 1. # magnetization oversaturated!!!
                print >>f, '  starting_magnetization(%d)=%s,' % (spcount,num2str(float(mag)))
                spcount += 1
        elif self.noncollinear:
            print >>f, '  noncolin=.true.,'
            if self.spinorbit:
                print >>f, '  lspinorb=.true.'
            spcount  = 1
            if self.nel == None:
                self.nvalence, self.nel = self.get_nvalence()
            for species in self.species: # FOLLOW SAME ORDERING ROUTINE AS FOR PSP                
                magindex = int(string.join([i for i in species if i.isdigit()],''))
                el  = species.strip('0123456789')
                mag = self.specdict[el]['magmoms'][magindex-1]/self.nel[el]
                assert np.abs(mag) <= 1. # magnetization oversaturated!!!
                print >>f, '  starting_magnetization(%d)=%s,' % (spcount,num2str(float(mag)))
                spcount += 1
        print >>f, '  input_dft=\''+self.xc+'\','
        edir = 3
        if dipfield:
            try:
                edir = self.dipole['edir']
            except:
                pass
        elif efield:
            try:
                edir = self.field['edir']
            except:
                pass
        if dipfield or efield:
            print >>f, '  edir='+str(edir)+','
        if dipfield:
            if self.dipole.has_key('emaxpos'):
                emaxpos = self.dipole['emaxpos']
            else:
                emaxpos = self.find_max_empty_space(edir)
            if self.dipole.has_key('eopreg'):
                eopreg = self.dipole['eopreg']
            else:
                eopreg = 0.025
            if self.dipole.has_key('eamp'):
                eamp = self.dipole['eamp']
            else:
                eamp = 0.0
            print >>f, '  emaxpos='+num2str(emaxpos)+','
            print >>f, '  eopreg='+num2str(eopreg)+','
            print >>f, '  eamp='+num2str(eamp)+','
        if efield:
            if self.field.has_key('emaxpos'):
                emaxpos = self.field['emaxpos']
            else:
                emaxpos = 0.0
            if self.field.has_key('eopreg'):
                eopreg = self.field['eopreg']
            else:
                eopreg = 0.0
            if self.field.has_key('eamp'):
                eamp = self.field['eamp']
            else:
                eamp = 0.0
            print >>f, '  emaxpos2='+num2str(emaxpos)+','
            print >>f, '  eopreg2='+num2str(eopreg)+','
            print >>f, '  eamp2='+num2str(eamp)+','
        if self.U is not None or self.J is not None or self.U_alpha is not None:
            print >>f, '  lda_plus_u=.true.,'
            if self.J is not None:
                print >>f, '  lda_plus_u_kind=1,'
            else:
                print >>f, '  lda_plus_u_kind=0,'
            print >>f, '  U_projection_type=\"%s\",' % (self.U_projection_type)
            if self.U is not None:
                for i,s in enumerate(self.species):
                    el = s.strip('0123456789')
                    if self.U.has_key(i):
                        Ui = self.U[i]
                    elif self.U.has_key(el):
                        Ui = self.U[el]
                    else:
                        Ui = 0.0
                    print >>f, '  Hubbard_U('+str(i+1)+')='+num2str(Ui)+','
            if self.J is not None:
                for i,s in enumerate(self.species):
                    el = s.strip('0123456789')
                    Ji ='KK'
                    if self.U.has_key(i):
                         Ji = self.J[i]
                    elif self.U.has_key(el):
                         Ji = self.J[el]
                    if Ji != 'KK':
                        print >>f, '  Hubbard_J('+str(i+1)+')='+num2str(Ji)+','
            if self.U_alpha is not None:
                for i, s in enumerate(self.species):
                    el = s.strip('0123456789')
                    if self.U_alpha.has_key(i):
                         U_alphai = self.U_alpha[i]
                    elif self.U.has_key(el):
                         U_alphai = self.U_alpha[el]
                    else:
                         U_alphai = 0.0
                    print >>f, '  Hubbard_alpha('+str(i+1)+')='+num2str(U_alphai)+','
        
        if self.nosym:
            print >>f, '  nosym=.true.,'
        if self.noinv:
            print >>f, '  noinv=.true.,'
        if self.nosym_evc:
            print >>f, '  nosym_evc=.true.,'
        if self.no_t_rev:
            print >>f, '  no_t_rev=.true.,'
        
        ### &ELECTRONS ###
        print >>f,'/\n&ELECTRONS'
        try:
            diag = self.convergence['diag']
            print >>f,'  diagonalization=\''+diag+'\','
        except:
            pass
        
        if self.calcmode!='hund':
            print >>f, '  conv_thr='+num2str(self.conv_thr)+','
        else:
            print >>f, '  conv_thr='+num2str(self.conv_thr*500.)+','
        for x in self.convergence.keys():
            if x=='mixing':
                print >>f, '  mixing_beta='+num2str(self.convergence[x])+','
            elif x=='maxsteps':
                print >>f, '  electron_maxstep='+str(self.convergence[x])+','
            elif x=='nmix':
                print >>f, '  mixing_ndim='+str(self.convergence[x])+','
            elif x=='mixing_mode':
                print >>f, '  mixing_mode=\''+self.convergence[x]+'\','
            elif x=='diago_cg_maxiter':
                print >>f, '  diago_cg_maxiter='+str(self.convergence[x])+','
        if self.startingpot is not None and self.calcmode!='hund':
            print >>f, '  startingpot=\''+self.startingpot+'\','
        if self.startingwfc is not None and self.calcmode!='hund':
            print >>f, '  startingwfc=\''+self.startingwfc+'\','

        ### &IONS ###
        if self.opt_algorithm is not None and self.calcmode not in ('scf','hund'):
            if self.cancalc:
                print >>f, '/\n&IONS\n  ion_dynamics=\''+self.opt_algorithm+'\','
            else:
                print >>f, '/\n&IONS\n  ion_dynamics=\'bfgs\','

        ### &CELL ###
        if self.cell_dynamics is not None:
            print >>f, '/\n&CELL\n  cell_dynamics=\''+self.cell_dynamics+'\','
            if self.press is not None:
                print >>f, '  press='+num2str(self.press)+','
            if self.dpress is not None:
                print >>f, '  press_conv_thr='+num2str(self.dpress)+','
            if self.cell_factor is not None:
                print >>f, '  cell_factor='+num2str(self.cell_factor)+','

        ### CELL_PARAMETERS
        print >>f, '/\nCELL_PARAMETERS'
        for i in range(3):
            print >>f, '%21.15fd0 %21.15fd0 %21.15fd0' % tuple(self.atoms.cell[i])

        print >>f, 'ATOMIC_SPECIES'
        for species in self.species:   # PSP ORDERING FOLLOWS SPECIESINDEX
            el = species.strip('0123456789')
            print >>f, species, self.specdict[el]['mass'][0], el+'.UPF'
        
        print >>f, 'ATOMIC_POSITIONS {crystal}'
        for species, mass, magmom, pos in self.specprops:
            print >>f, '%-4s %21.15fd0 %21.15fd0 %21.15fd0' % (species,pos[0],pos[1],pos[2])
        
        print >>f, 'K_POINTS automatic'
        if overridekpts is None:
            print >>f, self.kpts[0], self.kpts[1],self.kpts[2],
        else:
            print >>f, overridekpts[0], overridekpts[1],overridekpts[2],
        if overridekptshift is None:
            print >>f, self.kptshift[0],self.kptshift[1],self.kptshift[2]
        else:
            print >>f, overridekptshift[0],overridekptshift[1],overridekptshift[2]
        ### closing PWscf input file ###
        f.close()
        if self.verbose == 'high':
            print '\nPWscf input file %s written\n' % fname
        
    def set_atoms(self, atoms):
        if self.atoms is None or not self.started:
            self.atoms = atoms.copy()
        else:
            msg = 'stop calculator before assigning new atoms to it by using calc.stop()'
            if len(atoms)!=len(self.atoms):
                raise ValueError, msg
            
            x = atoms.cell-self.atoms.cell
            if np.max(x)>1E-13 or np.min(x)<-1E-13 or \
                (atoms.get_atomic_numbers()!=self.atoms.get_atomic_numbers()).any():
                raise ValueError, msg
            x = atoms.positions-self.atoms.positions
            if np.max(x)>1E-13 or np.min(x)<-1E-13 or (not self.started and not self.got_energy):
                self.recalculate = True
        self.atoms = atoms.copy()

    def update(self, atoms):
        if self.atoms is None:
            self.set_atoms(atoms)
        x = atoms.positions-self.atoms.positions
        if np.max(x)>1E-13 or np.min(x)<-1E-13 or (not self.started and not self.got_energy) or self.recalculate:
            self.recalculate = True
            if self.opt_algorithm!='ase3':
                self.stop()
            self.read(atoms)
        elif self.only_init:
            self.read(atoms)
        else:
            self.atoms = atoms.copy()

    def get_name(self):
        return 'QE-ASE3 interface'

    def get_version(self):
        return '0.1'

    def get_stress(self, atoms):
        raise NotImplementedError, 'stress interface not implemented\ntry using QE\'s internal relaxation routines instead'

    def init_only(self, atoms):
        if self.atoms is None:
            self.set_atoms(atoms)
        x = atoms.positions-self.atoms.positions
        if np.max(x)>1E-13 or np.min(x)<-1E-13 or (not self.started and not self.got_energy) or self.recalculate:
            self.recalculate = True
            if self.opt_algorithm!='ase3':
                self.stop()

            if not self.started:
                self.initialize(atoms)
                self.only_init = True
            elif self.recalculate:
                self.only_init = True
                if self.opt_algorithm == 'ase3':
                    p = atoms.positions
                    self.atoms = atoms.copy()
                    print >>self.cinp, 'G'
                    for x in p:
                        print >>self.cinp, ('%.15e %.15e %.15e' % (x[0],x[1],x[2])).replace('e','d')
                self.cinp.flush()

    def read(self, atoms):
        if self.writeversion:
            self.writeversion = False
            s = open(self.log,'a')
            s.write('  python dir         : '+self.mypath+'\n')
            exedir = os.path.dirname(os.popen('which pw.x').readline())
            s.write('  espresso dir       : '+exedir+'\n')
            s.write('  pseudo dir         : '+self.psppath+'\n')
            s.write('  espresso py svn    : '+svnver+'\n\n\n')
            s.close()

        if not self.started and not self.only_init:
            fresh = True
            self.initialize(atoms)
        else:
            fresh = False
        if self.recalculate:
            if not fresh and not self.only_init:
                if self.opt_algorithm == 'ase3':
                    p = atoms.positions
                    self.atoms = atoms.copy()
                    print >>self.cinp, 'G'
                    for x in p:
                        print >>self.cinp, ('%.15e %.15e %.15e' % (x[0],x[1],x[2])).replace('e','d')
                self.cinp.flush()
            self.only_init = False
            s = open(self.log,'a')
            a = self.cout.readline()
            s.write(a)
            atom_occ = {}
            while a!='' and a[:17]!='!    total energy' and a[:13]!='     stopping' and a[:20]!='     convergence NOT': 
                a = self.cout.readline()
                s.write(a)
                s.flush()
                
                if a[:19]=='     iteration #  1':
                    while (a!='' and a[:17]!='!    total energy' and a[:13]!='     stopping' and a[:20]!='     convergence NOT' and 
                           a[:22]!=' --- exit write_ns ---' ) :
                        a = self.cout.readline()   
                        s.write(a)
                        s.flush()
                        if a[:5]=='atom ':
                            atomnum = int(a[8:10])
                            if a[12:25]=='Tr[ns(na)] = ':#'atom    1   Tr[ns(na)] =   1.00000'
                                N0 = float(a[27:35])/2.
                            elif a[12:42]=='Tr[ns(na)] (up, down, total) =':
                                #'   4.20435  1.27943  5.48377'
                                N0 = [float(a[42:52]), float(a[53:62]), float(a[63:71])]
                                N0=N0[-1] # only taking the total occupation
                            atom_occ[atomnum-1]={}
                            atom_occ[atomnum-1][0]=N0
                if a[:39]=='     End of self-consistent calculation':
                    while a!='' and a[:17]!='!    total energy' and a[:13]!='     stopping' and a[:20]!='     convergence NOT':
                        a = self.cout.readline()
                        s.write(a) 
                        s.flush()
                        if a[:5]=='atom ':
                            atomnum = int(a[8:10])
                            if a[12:25]=='Tr[ns(na)] = ':#'atom    1   Tr[ns(na)] =   1.00000'
                                Nks = float(a[27:35])/2.
                            elif a[12:42]=='Tr[ns(na)] (up, down, total) =':
                                #'   4.20435  1.27943  5.48377'
                                Nks = [float(a[42:52]), float(a[53:62]), float(a[63:71])]
                                Nks=Nks[-1] # only taking the total occupation
                            atom_occ[atomnum-1]['ks']=Nks
                    break
            if a[:20]=='     convergence NOT':
                self.stop()
                raise RuntimeError, 'scf cycles did not converge\nincrease maximum number of steps and/or decreasing mixing'
            elif a[:13]=='     stopping':
                self.stop()
                self.checkerror()
                #if checkerror shouldn't find an error here,
                #throw this generic error
                raise RuntimeError, 'SCF calculation failed'
            elif a=='' and self.calcmode in ('ase3','relax','scf','vc-relax','vc-md','md'):
                self.checkerror()
                #if checkerror shouldn't find an error here,
                #throw this generic error
                raise RuntimeError, 'SCF calculation didn\'t converge'
            self.atom_occ = atom_occ
            if self.calcmode in ('ase3','relax','scf','vc-relax','vc-md','md','hund'):
                self.energy_free = float(a.split()[-2])*rydberg
                # get S*T correction (there is none for Marzari-Vanderbilt=Cold smearing)
                if self.occupations=='smearing' and self.calcmode!='hund' and self.smearing[0].upper()!='M' and self.smearing[0].upper()!='C':
                    a = self.cout.readline()
                    s.write(a)
                    while a[:13]!='     smearing':
                        a = self.cout.readline()
                        s.write(a)
                    self.ST = -float(a.split()[-2])*rydberg
                    self.energy_zero = self.energy_free + 0.5*self.ST
                else:
                    self.energy_zero = self.energy_free    
            else:
                self.energy_free = None
                self.energy_zero = None

            self.got_energy = True

            a = self.cout.readline()
            s.write(a)
            s.flush()

            if self.calcmode in ('ase3','relax','scf','vc-relax','vc-md','md'):
                if self.opt_algorithm == 'ase3' and self.calcmode != 'scf':
                    sys.stdout.flush()
                    while a[:5]!=' !ASE':
                        a = self.cout.readline()
                        s.write(a)
                        s.flush()
                    if not hasattr(self, 'forces'):
                        self.forces = np.empty((self.natoms,3), np.float)
                    for i in range(self.natoms):
                        self.cout.readline()
                    for i in range(self.natoms):
                        self.forces[i][:] = [float(x) for x in self.cout.readline().split()]
                    self.forces *= rydberg_over_bohr
                else:
                    a = self.cout.readline()
                    s.write(a)
                    while a[:11]!='     Forces':
                        a = self.cout.readline()
                        s.write(a)
                        s.flush()
                    a = self.cout.readline()
                    s.write(a)
                    if not hasattr(self, 'forces'):
                        self.forces = np.empty((self.natoms,3), np.float)
                    for i in range(self.natoms):
                        a = self.cout.readline()
                        while a.find('force')<0:
                            s.write(a) 
                            a = self.cout.readline()
                        s.write(a)
                        forceinp = a.split()
                        self.forces[i][:] = [float(x) for x in forceinp[len(forceinp)-3:]]
                    self.forces *= rydberg_over_bohr
            else:
                self.forces = None
            self.recalculate = False
            s.close()
            if self.opt_algorithm != 'ase3':
                self.stop()

            #get final energy and forces for internal QE relaxation run
            if self.calcmode in ('relax','vc-relax','vc-md','md'):
                if self.opt_algorithm == 'ase3':
                    self.stop()
                p = os.popen('grep -n "!    total" '+self.log+' | tail -1','r')
                n = int(p.readline().split(':')[0])-1
                p.close()
                f = open(self.log,'r')
                for i in range(n):
                    f.readline()
                self.energy_free = float(f.readline().split()[-2])*rydberg
                # get S*T correction (there is none for Marzari-Vanderbilt=Cold smearing)
                if self.occupations=='smearing' and self.calcmode!='hund' and self.smearing[0].upper()!='M' and self.smearing[0].upper()!='C':
                    a = f.readline()
                    while a[:13]!='     smearing':
                        a = f.readline()
                    self.ST = -float(a.split()[-2])*rydberg
                    self.energy_zero = self.energy_free + 0.5*self.ST
                else:
                    self.energy_zero = self.energy_free

                if self.U_projection_type == 'atomic':
                        a = f.readline()
                        while a[:11]!='     Forces':
                            a = f.readline()
                        f.readline()
                        if not hasattr(self, 'forces'):
                            self.forces = np.empty((self.natoms,3), np.float)
                        for i in range(self.natoms):
                            a = f.readline()
                            while a.find('force')<0:
                                a = f.readline()
                            forceinp = a.split()
                            self.forces[i][:] = [float(x) for x in forceinp[len(forceinp)-3:]]
                        self.forces *= rydberg_over_bohr                    
                f.close()
                
            self.checkerror()
    

    def initialize(self, atoms):
        """ Create the pw.inp input file and start the calculation. 
        If onlycreatepwinp is specified in calculator setup,
        only the input file will be written for manual submission.
        """
        if not self.started:
            self.atoms = atoms.copy()
            
            self.atoms2species()
            #s = a.get_chemical_symbols()
            #m = a.get_masses()
            #sd = {}
            #for x in zip(s, m):
            #    sd[x[0]] = x[1]
            #k = sd.keys()
            #k.sort()
            #self.species = [(x,sd[x]) for x in k] # UPDATE: NOT COMPATIBLE WITH MAGNETIC PARTS
            #self.nspec = len(self.species)
            self.natoms = len(self.atoms)
            #self.spos = zip(s, a.get_scaled_positions()) # UPDATE to have species indices
            self.check_spinpol()
            self.writeinputfile()
        if self.cancalc:
            self.start()

    def check_spinpol(self):
        mm = self.atoms.get_initial_magnetic_moments()
        sp = mm.any()
        self.summed_magmoms = np.sum(mm)
        if sp:
            if not self.spinpol and not self.noncollinear:
                raise KeyError('Explicitly specify spinpol=True or noncollinear=True for spin-polarized systems')
            elif abs(self.sigma) <= self.sigma_small and not self.fix_magmom:
                raise KeyError('Please use fix_magmom=True for sigma=0.0 eV and spinpol=True. Hopefully this is not an extended system...?')
        else:
            if self.spinpol and abs(self.sigma) <= self.sigma_small:
                self.fix_magmom = True
        if abs(self.sigma) <= self.sigma_small:
            self.occupations = 'fixed'

    def start(self):
        if not self.started:
	    if self.single_calculator:
		while len(espresso_calculators)>0:
		    espresso_calculators.pop().stop()
		espresso_calculators.append(self)
            if site.batch:
                cdir = os.getcwd()
                os.chdir(self.localtmp)
                os.system(site.perHostMpiExec+' cp '+self.localtmp+'/pw.inp '+self.scratch)
                if self.calcmode!='hund':
                    if not self.proclist:
                        self.cinp, self.cout = os.popen2(site.perProcMpiExec+' -wdir '+self.scratch+' pw.x '+self.parflags+' -in pw.inp')
                    else:
                        self.cinp, self.cout, self.cerr = os.popen3((site.perSpecProcMpiExec % (self.mycpus,self.myncpus))+' -wdir '+self.scratch+' pw.x '+self.parflags+' -in pw.inp|'+self.mypath+'/espfilter '+str(self.natoms)+' '+self.log+'0')
                else:
                    os.system(site.perProcMpiExec+' -wdir '+self.scratch+' pw.x -in pw.inp >>'+self.log)
                    os.system("sed s/occupations.*/occupations=\\'fixed\\',/ <"+self.localtmp+"/pw.inp | sed s/ELECTRONS/ELECTRONS\\\\n\ \ startingwfc=\\'file\\',\\\\n\ \ startingpot=\\'file\\',/ | sed s/conv_thr.*/conv_thr="+num2str(self.conv_thr)+",/ | sed s/tot_magnetization.*/tot_magnetization="+num2str(self.totmag)+",/ >"+self.localtmp+"/pw2.inp")
                    os.system(site.perHostMpiExec+' cp '+self.localtmp+'/pw2.inp '+self.scratch)
                    self.cinp, self.cout = os.popen2(site.perProcMpiExec+' -wdir '+self.scratch+' pw.x '+self.parflags+' -in pw2.inp')
                os.chdir(cdir)
            else:
                os.system('cp '+self.localtmp+'/pw.inp '+self.scratch)
                if self.calcmode!='hund':
                    self.cinp, self.cout = os.popen2('cd '+self.scratch+' ; '+'pw.x -in pw.inp')
                else:
                    os.system('cd '+self.scratch+' ; '+' pw.x -in pw.inp >>'+self.log)
                    os.system("sed s/occupations.*/occupations=\\'fixed\\',/ <"+self.localtmp+"/pw.inp | sed s/ELECTRONS/ELECTRONS\\\\n\ \ startingwfc=\\'file\\',\\\\n\ \ startingpot=\\'file\\',/ | sed s/conv_thr.*/conv_thr="+num2str(self.conv_thr)+",/ | sed s/tot_magnetization.*/tot_magnetization="+num2str(self.totmag)+",/ >"+self.localtmp+"/pw2.inp")
                    os.system('cp '+self.localtmp+'/pw2.inp '+self.scratch)
                    self.cinp, self.cout = os.popen2('cd '+self.scratch+' ; '+'pw.x -in pw2.inp')

            self.started = True

    def stop(self):
        if self.started:
            if self.opt_algorithm == 'ase3':
                #sending 'Q' to espresso tells it to quit cleanly
                print >>self.cinp, 'Q'
                try:
                    self.cinp.flush()
                except IOError:
                    #espresso may have already shut down, so flush may fail
                    pass
            else:
                self.cinp.flush()
            s = open(self.log,'a')
            a = self.cout.readline()
            s.write(a)
            while a!='':
                a = self.cout.readline()
                s.write(a)
                s.flush()
            s.close()
            self.cinp.close()
            self.cout.close()
            self.started = False


    def topath(self, filename):
        if os.path.isabs(filename):
            return filename
        else:
            return os.path.join(self.sdir,filename)


    def save_output(self, filename='calc.tgz'):
        file = self.topath(filename)
        self.update(self.atoms)
        self.stop()
        
        os.system('tar czf '+filename+' --directory='+self.scratch+' calc.save')


    def load_output(self, filename='calc.tgz'):
        self.stop()
        file = self.topath(filename)

        os.system('tar xzf '+filename+' --directory='+self.scratch)


    def save_wf(self, filename='wf.tgz'):
        file = self.topath(filename)
        self.update(self.atoms)
        self.stop()
        
        os.system('tar czf '+filename+' --directory='+self.scratch+' --exclude=calc.save .')


    def load_wf(self, filename='wf.tgz'):
        self.stop()
        file = self.topath(filename)

        os.system('tar xzf '+filename+' --directory='+self.scratch)


    def get_final_structure(self):
        """
        returns Atoms object according to a structure
        optimized internally by quantum espresso
        """
        from ase import Atoms

        self.stop()

        p = os.popen('grep -n Giannozzi '+self.log+'| tail -1','r')
        n = int(p.readline().split()[0].strip(':'))
        p.close()
        
        s = open(self.log,'r')
        #skip over previous runs in log in case the current log has been
        #appended to old ones
        for i in range(n):
            s.readline()
        
        a = s.readline()
        while a[:11]!='     celldm':
            a = s.readline()
        alat = float(a.split()[1])/1.889726
        a = s.readline()
        while a[:12]!='     crystal':
            a = s.readline()
        cell = []
        for i in range(3):
            cell.append([float(x) for x in s.readline().split()[3:6]])
        cell = np.array(cell)
        a = s.readline()
        while a[:12]!='     site n.':
            a = s.readline()
        pos = []
        syms = ''
        y = s.readline().split()
        while len(y)>0:
            nf = len(y)
            pos.append([float(x) for x in y[nf-4:nf-1]])
            syms += y[1].strip('0123456789')
            y = s.readline().split()
        pos = np.array(pos)*alat
        natoms = len(pos)
        
        #create atoms object with coordinates and unit cell
        #as specified in the initial ionic step in log
        atoms = Atoms(syms, pos, cell=cell*alat, pbc=(1,1,1))
        
        coord = 'angstrom)'
        a = s.readline()
        while a!='':
            while a[:7]!='CELL_PA' and a[:7]!='ATOMIC_' and a!='':
                a = s.readline()
            if a=='':
                break
            if a[0]=='A':
                coord = a.split('(')[-1]
                for i in range(natoms):
                    pos[i][:] = s.readline().split()[1:]
            else:
                for i in range(3):
                    cell[i][:] = s.readline().split()
            a = s.readline()

        atoms.set_cell(cell*alat, scale_atoms=False)
        
        if coord=='alat)':
            atoms.set_positions(pos*alat)
        elif coord=='bohr)':
            atoms.set_positions(pos*bohr)
        elif coord=='angstrom)':
            atoms.set_positions(pos)
        else:
            atoms.set_scaled_positions(pos)
        
        return atoms


    def get_final_stress(self):
        """
        returns 3x3 stress tensor after an internal
        unit cell relaxation in quantum espresso
        (also works for calcstress=True)
        """

        self.stop()
        
        p = os.popen('grep -3 "total   stress" '+self.log+' | tail -3','r')
        s = p.readlines()
        p.close()
        
        if len(s)!=3:
            raise RuntimeError, 'stress was not calculated\nconsider specifying calcstress or running a unit cell relaxation'
        
        stress = np.empty((3,3), np.float)
        for i in range(3):
            stress[i][:] = [float(x) for x in s[i].split()[:3]]
        
        return stress * rydberg/bohr**3


    def checkerror(self):
        p = os.popen('grep -n Giannozzi '+self.log+' | tail -1','r')
        try:
            n = int(p.readline().split()[0].strip(':'))
        except:
            raise RuntimeError, 'Espresso executable doesn\'t seem to have been started.'
        p.close()

        p = os.popen(('tail -n +%d ' % n)+self.log+' | grep -n %%%%%%%%%%%%%%%% |tail -2','r')
        s = p.readlines()
        p.close()

        if len(s)<2:
            return
        
        a = int(s[0].split()[0].strip(':'))+1
        b = int(s[1].split()[0].strip(':'))-a
        
        if b<1:
            return
        
        p = os.popen(('tail -n +%d ' % (a+n-1))+self.log+('|head -%d' % b),'r')
        err = p.readlines()
        p.close()
        
        if err[0].lower().find('error')<0:
            return
        
        msg = ''
        for e in err:
            msg += e
        raise RuntimeError, msg[:len(msg)-1]
   
    def relax_cell_and_atoms(self, 
            cell_dynamics='bfgs', # {'none', 'sd', 'damp-pr', 'damp-w', 'bfgs'}
            opt_algorithm='bfgs', # {'bfgs', 'damp'}
            cell_factor=1.2,
            fmax=None,
            press=None,
            dpress=None
            ):
        self.stop()
        oldmode = self.calcmode
        oldalgo = self.opt_algorithm
        oldcell = self.cell_dynamics
        oldfactor = self.cell_factor
        self.cell_dynamics=cell_dynamics
        self.opt_algorithm=opt_algorithm
        self.cell_factor=cell_factor
        oldfmax = self.fmax
        oldpress = self.press
        olddpress = self.dpress
        
        if fmax is not None:
            self.fmax = fmax
        if press is not None:
            self.press = press
        if dpress is not None:
            self.dpress = dpress
        self.calcmode='vc-relax'
        self.recalculate=True
        self.read(self.atoms)
        self.calcmode = oldmode
        self.opt_algorithm = oldalgo
        self.cell_dynamics = oldcell
        self.cell_factor = oldfactor
        self.fmax = oldfmax
        self.press = oldpress
        self.dpress = olddpress
   
    def relax_atoms(self,
            opt_algorithm='bfgs', # {'bfgs', 'damp'}
            fmax=None
            ):
        self.stop()
        oldmode = self.calcmode
        oldalgo = self.opt_algorithm
        self.opt_algorithm=opt_algorithm
        oldfmax = self.fmax
       
        self.calcmode='relax'
        if fmax is not None:
            self.fmax = fmax
        self.recalculate=True
        self.read(self.atoms)
        self.calcmode = oldmode
        self.opt_algorithm = oldalgo
        self.fmax = oldfmax


    #runs one of the .x binaries of the espresso suite
    #inp is expected to be in self.localtmp
    #log will be created in self.localtmp
    def run_espressox(self, binary, inp, log=None, piperead=False,
        parallel=True):
        if log is None:
            ll = ''
        else:
            ll = ' >>'+self.localtmp+'/'+log
        if site.batch and parallel:
            cdir = os.getcwd()
            os.chdir(self.localtmp)
            os.system(site.perHostMpiExec+' cp '+self.localtmp+'/'+inp+' '+self.scratch)
            if piperead:
                p = os.popen(site.perProcMpiExec+' -wdir '+self.scratch+' '+binary+' '+self.parflags+' -in '+inp+ll, 'r')
            else:
                os.system(site.perProcMpiExec+' -wdir '+self.scratch+' '+binary+' '+self.parflags+' -in '+inp+ll)
            os.chdir(cdir)
        else:
            os.system('cp '+self.localtmp+'/'+inp+' '+self.scratch)
            if piperead:
                p = os.popen('cd '+self.scratch+' ; '+binary+' -in '+inp+ll)
            else:
                os.system('cd '+self.scratch+' ; '+binary+' -in '+inp+ll)
        if piperead:
            return p

    def run_ppx(self, inp, log=None, inputpp=[], plot=[],
        output_format=5, iflag=3, piperead=False, parallel=True):
        if self.output.has_key('disk_io'):
            if self.output['disk_io'] == 'none':
                print "run_ppx requires output['disk_io'] to be at least 'low' and avoidio=False"
        self.stop()
        f = open(self.localtmp+'/'+inp, 'w')
        print >>f, '&INPUTPP\n  prefix=\'calc\',\n  outdir=\'.\','
        for a,b in inputpp:
            if type(b)==float:
                c = num2str(b)
            elif type(b)==str:
                c = "'"+b+"'"
            else:
                c = str(b)
            print >>f, '  '+a+'='+c+','
        print >>f, '/'
        print >>f, '&PLOT\n  iflag=%d,\n  output_format=%d,' % (iflag,output_format)
        for a,b in plot:
            if type(b)==float:
                c = num2str(b)
            elif type(b)==str:
                c = "'"+b+"'"
            else:
                c = str(b)
            print >>f, '  '+a+'='+c+','
        print >>f, '/'
        f.close()
        
        if piperead:
            return self.run_espressox('pp.x', inp, log=log,
                piperead=piperead, parallel=parallel)
        else:
            self.run_espressox('pp.x', inp, log=log, parallel=parallel)

    
    def get_fermi_level(self):
        self.stop()
        try:
            p = os.popen('grep Fermi '+self.log+'|tail -1', 'r')
            efermi = float(p.readline().split()[-2])
            p.close()
        except:
            raise RuntimeError, 'get_fermi_level called before DFT calculation was run'
        return efermi


    def calc_pdos(self,
        Emin = None,
        Emax = None,
        DeltaE = None,
        nscf = False,
        tetrahedra = False,
        slab = False,
        kpts = None,
        kptshift = None,
        ngauss = None,
        sigma = None,
        nscf_fermilevel=False,
        add_higher_channels=False):

        efermi = self.get_fermi_level()

        # run a nscf calculation with e.g. tetrahedra or more k-points etc.
        if nscf:
            self.writeinputfile(filename='pwnscf.inp',
                mode='nscf', usetetrahedra=tetrahedra, overridekpts=kpts,
                overridekptshift=kptshift, suppressforcecalc=True)
            self.run_espressox('pw.x', 'pwnscf.inp', 'pwnscf.log')
            if nscf_fermilevel:
                p = os.popen('grep Fermi '+self.localtmp+'/pwnscf.log|tail -1', 'r')
                efermi = float(p.readline().split()[-2])
                p.close()
        
        # remove old wave function projections
        os.system('rm -f '+self.scratch+'/*_wfc*')
        # create input for projwfc.x
        f = open(self.localtmp+'/pdos.inp', 'w')
        print >>f, '&PROJWFC\n  prefix=\'calc\',\n  outdir=\'.\','
        if Emin is not None:
            print >>f, '  Emin = '+num2str(Emin+efermi)+','
        if Emax is not None:
            print >>f, '  Emax = '+num2str(Emax+efermi)+','
        if DeltaE is not None:
            print >>f, '  DeltaE = '+num2str(DeltaE)+','
        if slab:
            print >>f, '  lslab = .true.,'
        if ngauss is not None:
            print >>f, '  ngauss = '+str(ngauss)+','
        if sigma is not None:
            print >>f, '  degauss = '+num2str(sigma/rydberg)+','
        print >>f, '/'
        f.close()
        # run projwfc.x
        self.run_espressox('projwfc.x', 'pdos.inp', 'pdos.log')
        
        # read in total density of states
        dos = np.loadtxt(self.scratch+'/calc.pdos_tot')
        if len(dos[0])>3:
            nspin = 2
        else:
            nspin = 1
        self.dos_energies = dos[:,0] - efermi
        self.dos_total = dos[:,1]
        npoints = len(self.dos_energies)
        
        channels = {'s':0, 'p':1, 'd':2, 'f': 3}
        # read in projections onto atomic orbitals
        self.pdos = [{} for i in range(self.natoms)]
        p = os.popen('ls '+self.scratch+'/calc.pdos_atm*')
        proj = p.readlines()
        p.close()
        proj.sort()
        for i,inp in enumerate(proj):
            inpfile = inp.strip()
            pdosinp = np.genfromtxt(inpfile)
            spl = inpfile.split('#')
            iatom = int(spl[1].split('(')[0])-1
            channel = spl[2].split('(')[1].rstrip(')')
            #ncomponents = 2*l+1 +1  (latter for m summed up)
            ncomponents = (2*channels[channel]+2) * nspin
            if not self.pdos[iatom].has_key(channel):
                self.pdos[iatom][channel] = np.zeros((ncomponents,npoints), np.float)
                first = True
            else:
                first = False
            if add_higher_channels or first:
                for j in range(ncomponents):
                    self.pdos[iatom][channel][j] += pdosinp[:,(j+1)]
        
        return self.dos_energies, self.dos_total, self.pdos


    def read_3d_grid(self, stream, log):
        f = open(self.localtmp+'/'+log, 'a')
        x = stream.readline()
        while x!='' and x[:11]!='DATAGRID_3D':
            f.write(x)
            x = stream.readline()
        if x=='':
            raise RuntimeError, 'error reading 3D data grid'
        f.write(x)
        nx, ny, nz = [int(y) for y in stream.readline().split()]
        origin = np.array([float(y) for y in stream.readline().split()])
        cell = np.empty((3,3), np.float)
        for i in range(3):
            cell[i][:] = [float(y) for y in stream.readline().split()]
        data = np.reshape(np.fromfile(stream, count=nx*ny*nz, sep=' '),
            (nx,ny,nz), order='F')

        x = stream.readline()
        while x!='':
            f.write(x)
            x = stream.readline()

        f.close()
        return (origin,cell,data)


    def read_2d_grid(self, stream, log):
        f = open(self.localtmp+'/'+log, 'a')
        x = stream.readline()
        while x!='' and x[:11]!='DATAGRID_2D':
            f.write(x)
            x = stream.readline()
        if x=='':
            raise RuntimeError, 'error reading 2D data grid'
        f.write(x)
        nx, ny = [int(y) for y in stream.readline().split()]
        origin = np.array([float(y) for y in stream.readline().split()])
        cell = np.empty((3,3), np.float)
        for i in range(3):
            cell[i][:] = [float(y) for y in stream.readline().split()]
        data = np.reshape(np.fromfile(stream, count=nx*ny, sep=' '),
            (nx,ny), order='F')

        x = stream.readline()
        while x!='':
            f.write(x)
            x = stream.readline()

        f.close()
        return (origin,cell,data)


    def extract_charge_density(self, spin='both'):
        if spin=='both' or spin==0:
            s = 0
        elif spin=='up' or spin==1:
            s = 1
        elif spin=='down' or spin==2:
            s = 2
        else:
            raise ValueError, 'unknown spin component'
        
        p = self.run_ppx('charge.inp',
            inputpp=[['plot_num',0],['spin_component',s]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'charge.log')
        p.close()
        return (origin,cell,data)

    def xsf_charge_density(self, xsf, spin='both'):
        if spin=='both' or spin==0:
            s = 0
        elif spin=='up' or spin==1:
            s = 1
        elif spin=='down' or spin==2:
            s = 2
        else:
            raise ValueError, 'unknown spin component'
        
        self.run_ppx('charge.inp',
            inputpp=[['plot_num',0],['spin_component',s]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='charge.log')


    def extract_total_potential(self, spin='both'):
        if spin=='both' or spin==0:
            s = 0
        elif spin=='up' or spin==1:
            s = 1
        elif spin=='down' or spin==2:
            s = 2
        else:
            raise ValueError, 'unknown spin component'
        
        p = self.run_ppx('totalpot.inp',
            inputpp=[['plot_num',1],['spin_component',s]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'totalpot.log')
        p.close()
        return (origin,cell,data)

    def xsf_total_potential(self, xsf, spin='both'):
        if spin=='both' or spin==0:
            s = 0
        elif spin=='up' or spin==1:
            s = 1
        elif spin=='down' or spin==2:
            s = 2
        else:
            raise ValueError, 'unknown spin component'
        
        self.run_ppx('totalpot.inp',
            inputpp=[['plot_num',1],['spin_component',s]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='totalpot.log')


    def extract_local_ionic_potential(self):
        p = self.run_ppx('vbare.inp',
            inputpp=[['plot_num',2]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'vbare.log')
        p.close()
        return (origin,cell,data)

    def xsf_local_ionic_potential(self, xsf):
        self.run_ppx('vbare.inp',
            inputpp=[['plot_num',2]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='vbare.log')


    def extract_local_dos_at_efermi(self):
        p = self.run_ppx('ldosef.inp',
            inputpp=[['plot_num',3]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'ldosef.log')
        p.close()
        return (origin,cell,data)

    def xsf_local_dos_at_efermi(self, xsf):
        self.run_ppx('ldosef.inp',
            inputpp=[['plot_num',3]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='ldosef.log')


    def extract_local_entropy_density(self):
        p = self.run_ppx('lentr.inp',
            inputpp=[['plot_num',4]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'lentr.log')
        p.close()
        return (origin,cell,data)

    def xsf_local_entropy_density(self, xsf):
        self.run_ppx('lentr.inp',
            inputpp=[['plot_num',4]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='lentr.log')


    def extract_stm_data(self, bias):
        p = self.run_ppx('stm.inp',
            inputpp=[['plot_num',5],['sample_bias',bias/rydberg]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'stm.log')
        p.close()
        return (origin,cell,data)

    def xsf_stm_data(self, xsf, bias):
        self.run_ppx('stm.inp',
            inputpp=[['plot_num',5],['sample_bias',bias/rydberg]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='stm.log')


    def extract_magnetization_density(self):
        p = self.run_ppx('magdens.inp',
            inputpp=[['plot_num',6]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'magdens.log')
        p.close()
        return (origin,cell,data)

    def xsf_magnetization_density(self, xsf):
        self.run_ppx('magdens.inp',
            inputpp=[['plot_num',6]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='magdens.log')


    def extract_wavefunction_density(self, band, kpoint=0, spin='up',
        gamma_with_sign=False):
        if spin=='up' or spin==1:
            s = 0
        elif spin=='down' or spin==2:
            s = 1
        elif spin=='charge' or spin==0:
            s = 0
        elif spin=='x':
            s = 1
        elif spin=='y':
            s = 2
        elif spin=='z':
            s = 3
        else:
            raise ValueError, 'unknown spin component'
        if self.spinpol:
            p = os.popen('grep "number of k points=" '+self.log+'|tail -1|tr \'=\' \' \'', 'r')
            nkp = int(p.readline().split()[4])
            p.close()
            kp = kpoint+nkp/2*s
        else:
            kp = kpoint
        inputpp = [['plot_num',7],['kpoint',kp],['kband',band]]
        if gamma_with_sign:
            inputpp.append(['lsign','.true.'])
        if self.noncollinear:
            inputpp.append(['spin_component',s])
        p = self.run_ppx('wfdens.inp',
            inputpp=inputpp,
            piperead=True, parallel=True)
        origin,cell,data = self.read_3d_grid(p, 'wfdens.log')
        p.close()
        return (origin,cell,data)

    def xsf_wavefunction_density(self, xsf, band, kpoint=0, spin='up',
        gamma_with_sign=False):
        if spin=='up' or spin==1:
            s = 0
        elif spin=='down' or spin==2:
            s = 1
        elif spin=='charge' or spin==0:
            s = 0
        elif spin=='x':
            s = 1
        elif spin=='y':
            s = 2
        elif spin=='z':
            s = 3
        else:
            raise ValueError, 'unknown spin component'
        if self.spinpol:
            p = os.popen('grep "number of k points=" '+self.log+'|tail -1|tr \'=\' \' \'', 'r')
            nkp = int(p.readline().split()[4])
            p.close()
            kp = kpoint+nkp/2*s
        else:
            kp = kpoint
        inputpp = [['plot_num',7],['kpoint',kp],['kband',band]]
        if gamma_with_sign:
            inputpp.append(['lsign','.true.'])
        if self.noncollinear:
            inputpp.append(['spin_component',s])
        self.run_ppx('wfdens.inp',
            inputpp=inputpp,
            plot=[['fileout',self.topath(xsf)]],
            parallel=True, log='wfdens.log')


    def extract_electron_localization_function(self):
        p = self.run_ppx('elf.inp',
            inputpp=[['plot_num',8]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'elf.log')
        p.close()
        return (origin,cell,data)

    def xsf_electron_localization_function(self, xsf):
        self.run_ppx('elf.inp',
            inputpp=[['plot_num',8]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='elf.log')


    def extract_density_minus_atomic(self):
        p = self.run_ppx('dens_wo_atm.inp',
            inputpp=[['plot_num',9]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'dens_wo_atm.log')
        p.close()
        return (origin,cell,data)

    def xsf_density_minus_atomic(self, xsf):
        self.run_ppx('dens_wo_atm.inp',
            inputpp=[['plot_num',9]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='dens_wo_atm.log')


    def extract_int_local_dos(self, spin='both', emin=None, emax=None):
        if spin=='both' or spin==0:
            s = 0
        elif spin=='up' or spin==1:
            s = 1
        elif spin=='down' or spin==2:
            s = 2
        else:
            raise ValueError, 'unknown spin component'
        
        inputpp=[['plot_num',10],['spin_component',s]]
        efermi = self.get_fermi_level()
        if emin is not None:
            inputpp.append(['emin',emin-efermi])
        if emax is not None:
            inputpp.append(['emax',emax-efermi])
        
        p = self.run_ppx('ildos.inp',
            inputpp=inputpp,
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'ildos.log')
        p.close()
        return (origin,cell,data)

    def xsf_int_local_dos(self, xsf, spin='both', emin=None, emax=None):
        if spin=='both' or spin==0:
            s = 0
        elif spin=='up' or spin==1:
            s = 1
        elif spin=='down' or spin==2:
            s = 2
        else:
            raise ValueError, 'unknown spin component'
        
        inputpp=[['plot_num',10],['spin_component',s]]
        efermi = self.get_fermi_level()
        if emin is not None:
            inputpp.append(['emin',emin-efermi])
        if emax is not None:
            inputpp.append(['emax',emax-efermi])
        
        self.run_ppx('ildos.inp',
            inputpp=inputpp,
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='ildos.log')


    def extract_ionic_and_hartree_potential(self):
        p = self.run_ppx('potih.inp',
            inputpp=[['plot_num',11]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'potih.log')
        p.close()
        return (origin,cell,data)

    def xsf_ionic_and_hartree_potential(self, xsf):
        self.run_ppx('potih.inp',
            inputpp=[['plot_num',11]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='potih.log')


    def extract_sawtooth_potential(self):
        p = self.run_ppx('sawtooth.inp',
            inputpp=[['plot_num',12]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'sawtooth.log')
        p.close()
        return (origin,cell,data)

    def xsf_sawtooth_potential(self, xsf):
        self.run_ppx('sawtooth.inp',
            inputpp=[['plot_num',12]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='sawtooth.log')


    def extract_noncollinear_magnetization(self, spin='all'):
        if spin=='all' or spin=='charge' or spin==0:
            s = 0
        elif spin=='x':
            s = 1
        elif spin=='y':
            s = 2
        elif spin=='z':
            s = 3
        else:
            raise ValueError, 'unknown spin component'
        p = self.run_ppx('noncollmag.inp',
            inputpp=[['plot_num',13],['spin_component',s]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'noncollmag.log')
        p.close()
        return (origin,cell,data)

    def xsf_noncollinear_magnetization(self, xsf, spin='all'):
        if spin=='all' or spin=='charge' or spin==0:
            s = 0
        elif spin=='x':
            s = 1
        elif spin=='y':
            s = 2
        elif spin=='z':
            s = 3
        else:
            raise ValueError, 'unknown spin component'
        self.run_ppx('noncollmag.inp',
            inputpp=[['plot_num',13],['spin_component',s]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False)


    def extract_ae_charge_density(self, spin='both'):
        if spin=='both' or spin==0:
            s = 0
        elif spin=='up' or spin==1:
            s = 1
        elif spin=='down' or spin==2:
            s = 2
        else:
            raise ValueError, 'unknown spin component'
        
        p = self.run_ppx('aecharge.inp',
            inputpp=[['plot_num',17],['spin_component',s]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'aecharge.log')
        p.close()
        return (origin,cell,data)

    def xsf_ae_charge_density(self, xsf, spin='both'):
        if spin=='both' or spin==0:
            s = 0
        elif spin=='up' or spin==1:
            s = 1
        elif spin=='down' or spin==2:
            s = 2
        else:
            raise ValueError, 'unknown spin component'
        
        self.run_ppx('aecharge.inp',
            inputpp=[['plot_num',17],['spin_component',s]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='aecharge.log')


    def extract_noncollinear_xcmag(self):
        p = self.run_ppx('ncxcmag.inp',
            inputpp=[['plot_num',18]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'ncxcmag.log')
        p.close()
        return (origin,cell,data)

    def xsf_noncollinear_xcmag(self, xsf):
        self.run_ppx('ncxcmag.inp',
            inputpp=[['plot_num',18]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='ncxcmag.log')


    def extract_reduced_density_gradient(self):
        p = self.run_ppx('redgrad.inp',
            inputpp=[['plot_num',19]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'redgrad.log')
        p.close()
        return (origin,cell,data)

    def xsf_reduced_density_gradient(self, xsf):
        self.run_ppx('redgrad.inp',
            inputpp=[['plot_num',19]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='redgrad.log')


    def extract_middle_density_hessian_eig(self):
        p = self.run_ppx('mideig.inp',
            inputpp=[['plot_num',20]],
            piperead=True, parallel=False)
        origin,cell,data = self.read_3d_grid(p, 'mideig.log')
        p.close()
        return (origin,cell,data)

    def xsf_middle_density_hessian_eig(self, xsf):
        self.run_ppx('mideig.inp',
            inputpp=[['plot_num',20]],
            plot=[['fileout',self.topath(xsf)]],
            parallel=False, log='mideig.log')


    def find_max_empty_space(self, edir=3):
        """
        Assuming periodic boundary conditions, finds the largest
        continuous segment of free, unoccupied space and returns
        its midpoint in scaled coordinates (0 to 1) in the edir direction (default z).
        """
        position_array = self.atoms.get_scaled_positions()[..., edir - 1]  # 0-indexed direction
        position_array.sort()
        differences = np.diff(position_array)
        differences = np.append(differences, position_array[0] + 1 - position_array[-1])  # through the PBC
        max_diff_index = np.argmax(differences)
        if max_diff_index == len(position_array) - 1:
            return (position_array[0] + 1 + position_array[-1]) / 2.
        else:
            return (position_array[max_diff_index] + position_array[max_diff_index + 1]) / 2.


    def get_work_function(self, pot_filename="pot.xsf", edir=3):
        """
        Calculates the work function of a calculation by subtracting the electrostatic
        potential of the vacuum (from averaging the output of pp.x num_plot 11 in the z
        direction by default) from the Fermi energy.
        Values used for average.x come from the espresso example for work function for a surface
        TODO: Implement some sort of tuning for these parameters?
        """
        if pot_filename[0] != '/':
            file = self.sdir + '/' + pot_filename
        else:
            file = pot_filename
        self.update(self.atoms)
        self.stop()
        if not os.path.exists(file):
            self.run_ppx('wf_pp.in', log='wf_pp.log',
                inputpp=[('plot_num', 11), ('filplot', self.topath('pot.xsf'))], 
                output_format=3, iflag=3, piperead=False, parallel=False)

        f = open(self.localtmp + '/avg.in', 'w')
        print >>f, '1'
        print >>f, self.sdir + "/" + pot_filename
        print >>f, '1.D0'
        print >>f, '1440'
        print >>f, '3'
        print >>f, '3.835000000'
        print >>f, ''
        f.close()
        os.system('cp ' + self.localtmp + '/avg.in ' + self.scratch)
        os.system('cd ' + self.scratch + ' ; ' + 'average.x < avg.in >>' + self.localtmp + '/avg.out')

        # Pick a good place to sample vacuum level
        cell_length = self.atoms.cell[edir - 1][edir - 1] / bohr
        vacuum_pos = self.find_max_empty_space(edir) * cell_length
        avg_out = open(self.localtmp + '/avg.out', 'r')
        record = False
        average_data = []
        lines = list(avg_out)
        for line in lines:
            if len(line.split()) == 3 and line.split()[0] == "0.000000000":
                record = True
            elif len(line.split()) == 0:
                record = False
            if record == True:
                average_data.append([float(i) for i in line.split()])
        vacuum_energy = average_data[np.abs(np.array(average_data)[..., 0] - vacuum_pos).argmin()][2]

        # Get the latest Fermi energy
        fermi_data = os.popen('grep -n "Fermi" ' + self.log + ' | tail -1', 'r')
        fermi_energy = float(fermi_data.readline().split()[-2])
        fermi_data.close()

        return vacuum_energy * rydberg - fermi_energy


    def get_world(self):
        from worldstub import world
        return world(site.nprocs)
