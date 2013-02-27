#cluster-dependent definitions
scratch = '/scratch'
submitdir = '$LS_SUBCWD'
jobid = '$LSB_BATCH_JID'
getprocs = ' echo -e $LSB_HOSTS | sed s/" "/"\\\\n"/g >machinefile ;'\
          +' uniq machinefile >uniqmachinefile ;'\
          +' nodes=`wc -l <uniqmachinefile` ;'\
          +' np=`wc -l <machinefile` '

import os
fin,fout = os.popen4('mpirun --version')
mpiversion=fout.readlines()[0]
mpiversion=mpiversion.split()[3]
rsh_agent='orte_rsh_agent'
if mpiversion[:3]=='1.4':
        rsh_agent='plm_rsh_agent'
        
perHostMpiExec = 'mpiexec --mca '+rsh_agent+' /afs/slac.stanford.edu/package/lsf/bin.slac/gmmpirun_lsgrun.sh -machinefile uniqmachinefile -np `wc -l <uniqmachinefile`'
perProcMpiExec = 'pam -g /afs/slac/g/suncat/bin/suncat-tsmpirun -x LD_LIBRARY_PATH'

from ase.calculators.general import Calculator
import atexit
import sys, string
import numpy as np
from types import FileType, StringType


def checkbatch():
    p = os.popen('echo '+jobid, 'r')
    batch = (p.readline().strip()!='')
    p.close()
    return batch

def mklocaltmp(batch, odir):
    if batch:
        p = os.popen('echo '+submitdir,'r')
        s = p.readline().strip()
        p.close()
        job = jobid
    else:
        s = os.getcwd()
        job = ''
    if odir is None or len(odir)==0:
        p = os.popen('mktemp -d '+s+'/qe"'+job+'"_XXXXX', 'r')
        tdir = p.readline().strip()
        p.close()
    else:
        if odir[0]=='/':
            tdir = odir
        else:
            tdir = s+'/'+odir
        os.system('mkdir -p '+tdir)
    return tdir

def mkscratch(batch,localtmp):
    if batch:
        pernodeexec = perHostMpiExec
        job = jobid
    else:
        pernodeexec = ''
        job = ''
    p = os.popen('mktemp -d '+scratch+'/qe"'+job+'"_XXXXX', 'r')
    tdir = p.readline().strip()
    p.close()
    if pernodeexec!='':
        cdir = os.getcwd()
        os.chdir(localtmp)
        os.system(pernodeexec + ' mkdir -p '+tdir)
        os.chdir(cdir)
    return tdir

def mpisetup(tdir):
    cdir = os.getcwd()
    os.chdir(tdir)
    p = os.popen(getprocs+' ; sh -c "echo $nodes $np"', 'r')
    os.chdir(cdir)
    nodes,np = p.readline().split()
    p.close()
    return nodes,np

def cleanup(tmp, scratch, removewf, batch, calc):
    try:
        calc.stop()
    except:
        pass
    if batch:
        pernodeexec = perHostMpiExec
    else:
        pernodeexec = ''
    if removewf:
        os.system('rm -r '+scratch+'/*wfc* 2>/dev/null')
    os.system('cp -r '+scratch+' '+tmp)
    cdir = os.getcwd()
    os.chdir(tmp)
    os.system(pernodeexec + ' rm -r '+scratch+' 2>/dev/null')
    os.chdir(cdir)

def getsubmitorcurrentdir():
    p = os.popen('echo '+submitdir, 'r')
    s = p.readline().strip()
    p.close()
    if s!='':
        return s
    else:
        return os.getcwd()

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
        list1 == u_l1
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


hartree = 27.21138505
rydberg = 0.5*hartree
bohr = 0.52917721092
rydberg_over_bohr = rydberg / bohr

#add 'd0' to floating point number to avoid random trailing digits in Fortran input routines
def num2str(x):
    s = str(x)
    if s.find('e')<0:
        s += 'd0'
    return s


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
                 outdir = None,
                 calcstress = False,
                 smearing = 'fd',
                 sigma = 0.1,
                 fix_magmom = False,
                 U = None,
                 J = None,
                 U_alpha = None,
                 U_projection_type = 'atomic',
                 tot_charge = 0.0, # +1 means 1 e missing, -1 means 1 extra e
                 tot_magnetization = -1, #-1 means unspecified, 'hund' means Hund's rule for each atom
                 occupations = 'smearing', # 'smearing', 'fixed', 'tetrahedra'
                 dipole = {'status':False},
                 field = {'status':False},
                 output = {'avoidio':False,
                           'removewf':True,
                           'wf_collect':False},
                 convergence = {'energy':1e-6,
                                'mixing':0.7,
                                'maxsteps':100,
                                'diag':'david'},
                 startingpot = None,
                 startingwfc = None,
                 onlycreatepwinp = None, #specify filename to only create pw input
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
        self.atoms = None
        self.sigma_small = 1e-13
        self.started = False
        self.got_energy = False

        self.input_update() # Create the tmp output folder 


    def input_update(self):
        """ Run initialization functions, such that this can be called if variables in espresso are
        changes using set or directly. 
        """
        self.create_outdir() # Create the tmp output folder 

        #sdir is the directory the script is run or submitted from
        self.sdir = getsubmitorcurrentdir()

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


    def create_outdir(self):
        if self.onlycreatepwinp is None:
            self.batch = checkbatch()
            self.localtmp = mklocaltmp(self.batch, self.outdir)
            if self.batch:
                self.nodes,self.np = mpisetup(self.localtmp)
            self.scratch = mkscratch(self.batch, self.localtmp)
            if self.output is not None and self.output.has_key('removewf'):
                removewf = self.output['removewf']
            else:
                removewf = True
            atexit.register(cleanup, self.localtmp, self.scratch, removewf, self.batch, self)
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
                self.outdir = value
                self.create_outdir()
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
                self.pw = value
            if key == 'dw':
                self.dw = value
            if key == 'output':
                self.output = value
            if key == 'convergence':
                self.convergence = value
            if key == 'kpts':
                self.kpts = value
            if key == 'kshift':
                self.kshift = value
        self.input_update()

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
                            Ulist.append(0.)
                    if self.J is not None:
                        if i in self.J:
                            Jlist.append(self.J[i])
                        elif atom.symbol in self.J:                          
                            Jlist.append(self.J[atom.symbol])
                        else:
                            Jlist.append(0.)
                    if self.U_alpha is not None:
                        if i in self.U_alpha:
                            U_alphalist.append(self.U_alpha[i])
                        elif atom.symbol in self.U_alpha:
                            U_alphalist.append(self.U_alpha[atom.symbol])
                        else:
                            U_alphalist.append(0.)
            
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
           
    
    def writeinputfile(self):
        if self.atoms is None:
            raise ValueError, 'no atoms defined'
        if self.cancalc:
            fname = self.localtmp+'/pw.inp'
            #f = open(self.localtmp+'/pw.inp', 'w')
        else:
            fname = self.pwinp
            #f = open(self.pwinp, 'w')
        f = open(fname, 'w')
        
        ### &CONTROL ###
        if self.calcmode=='ase3':
            print >>f, '&CONTROL\n  calculation=\'relax\',\n  prefix=\'calc\','
        elif self.calcmode=='hund':
            print >>f, '&CONTROL\n  calculation=\'scf\',\n  prefix=\'calc\','
        else:
            print >>f, '&CONTROL\n  calculation=\''+self.calcmode+'\',\n  prefix=\'calc\','
        print >>f, '  pseudo_dir=\''+self.psppath+'\','
        print >>f, '  outdir=\'.\','
        efield = (self.field['status']==True)
        dipfield = (self.dipole['status']==True)
        if efield or dipfield:
            print >>f, '  tefield=.true.,'
            if dipfield:
                print >>f, '  dipfield=.true.,'
        if not self.dontcalcforces:
            print >>f, '  tprnfor=.true.,'
            if self.calcstress:
                print >>f, '  tstress=.true.,'
            if self.output is not None:
                if self.output.has_key('avoidio'):
                    if self.output['avoidio']:
                        print >>f, '  disk_io=\'none\','
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
        if not self.tot_charge:
            print >>f, '  tot_charge='+str(self.tot_charge)+','
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
                print >>f, '  nbnd='+str(self.nbands)+','
            else:
                n = 0
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
                for x in self.specprops:
                    n += nel[x[0].strip('0123456789')]
                if not self.spinpol:
                    n /= 2
                print >>f, '  nbnd='+str(n-self.nbands)+','
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
            for species in self.species: # FOLLOW SAME ORDERING ROUTINE AS FOR PSP                
                magindex = int(string.join([i for i in species if i.isdigit()],''))
                el  = species.strip('0123456789')
                mag = self.specdict[el]['magmoms'][magindex-1]
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
                emaxpos = 0.05
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
                        Ui = '1D-40'
                    if 'D' in Ui or 'd' in Ui or 'E' in Ui or 'e' in Ui:
                        print >>f, '  Hubbard_U('+str(i+1)+')='+str(Ui)+','
                    else:
                        print >>f, '  Hubbard_U('+str(i+1)+')='+str(Ui)+'d0,'
            if self.J is not None:
                for i,s in enumerate(self.species):
                    el = s.strip('0123456789')
                    Ji ='KK'
                    if self.U.has_key(i):
                         Ji = self.J[i]
                    elif self.U.has_key(el):
                         Ji = self.J[el]
                    if Ji != 'KK':
                        if 'D' in Ji or 'd' in Ji or 'E' in Ji or 'e' in Ji:
                            print >>f, '  Hubbard_J('+str(i+1)+')='+str(Ji)+','
                        else:
                            print >>f, '  Hubbard_J('+str(i+1)+')='+str(Ji)+'d0,'
            if self.U_alpha is not None:
                for i, s in enumerate(self.species):
                    el = s.strip('0123456789')
                    if self.U_alpha.has_key(i):
                         U_alphai = self.U_alpha[i]
                    elif self.U.has_key(el):
                         U_alphai = self.U_alpha[el]
                    else:
                         U_alphai = '1D-40'
                    if isinstance(U_alphai, str) and ('D' in U_alphai or 'd' in U_alphai or 
                        'E' in U_alphai or 'e' in U_alphai):
                         print >>f, '  Hubbard_alpha('+str(i+1)+')='+str(U_alphai)+','
                    else:
                         print >>f, '  Hubbard_alpha('+str(i+1)+')='+str(U_alphai)+'d0,'


                    """
                    if self.U_alpha.has_key(i):
                        print >>f, '  Hubbard_alpha('+str(i+1)+')='+str(self.U_alpha[i])+'d0,'
                    elif self.U_alpha.has_key(el):
                        print >>f, '  Hubbard_alpha('+str(i+1)+')='+str(self.U_alpha[el])+'d0,'
                    else:
                        print >>f, '  Hubbard_alpha('+str(i+1)+')=1D-40'
                    """
        
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
        if self.startingpot is not None and self.calcmode!='hund':
            print >>f, '  startingpot=\''+self.startingpot+'\','
        if self.startingwfc is not None and self.calcmode!='hund':
            print >>f, '  startingwfc=\''+self.startingwfc+'\','

        ### &IONS ###
        if self.opt_algorithm is not None:
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
        print >>f, self.kpts[0], self.kpts[1],self.kpts[2],self.kptshift[0],self.kptshift[1],self.kptshift[2]
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
            if max(x.flat)>1E-13 or min(x.flat)<-1E-13 or \
                atoms.get_atomic_numbers()!=self.atoms.get_atomic_numbers():
                raise ValueError, msg
        self.atoms = atoms.copy()

    def update(self, atoms):
        if self.atoms is None:
            self.set_atoms(atoms)
        x = atoms.positions-self.atoms.positions
        if max(x.flat)>1E-13 or min(x.flat)<-1E-13 or (not self.started and not self.got_energy):
            self.recalculate = True
            if self.opt_algorithm!='ase3':
                self.stop()
            self.read(atoms)
        else:
            self.atoms = atoms.copy()

    def get_name(self):
        return 'QE-ASE3 interface'

    def get_version(self):
        return '0.1'

    def get_stress(self, atoms):
        raise NotImplementedError, 'stress interface not implemented\ntry using QE\'s internal relaxation routines instead'

    def read(self, atoms):
        if not self.started:
            fresh = True
            self.initialize(atoms)
        else:
            fresh = False
        if self.recalculate:
            if not fresh:
                if self.opt_algorithm == 'ase3':
                    p = atoms.positions
                    self.atoms = atoms.copy()
                    print >>self.cinp, 'G'
                    for x in p:
                        print >>self.cinp, ('%.15e %.15e %.15e' % (x[0],x[1],x[2])).replace('e','d')
                self.cinp.flush()
            s = open(self.localtmp+'/log','a')
            a = self.cout.readline()
            s.write(a)
            atom_occ = {}
            while a!='' and a[:17]!='!    total energy' and a[:13]!='     stopping': 
                a = self.cout.readline()
                s.write(a)
                s.flush()
                
                if a[:19]=='     iteration #  1':
                    while (a!='' and a[:17]!='!    total energy' and a[:13]!='     stopping' and 
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
                    while a!='' and a[:17]!='!    total energy' and a[:13]!='     stopping':
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
            if a[:13]=='     stopping':
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
                p = os.popen('grep -n "!    total" '+self.localtmp+'/log | tail -1','r')
                n = int(p.readline().split(':')[0])-1
                p.close()
                f = open(self.localtmp+'/log','r')
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
            if not self.spinpol:
                raise KeyError('Explicitly specify spinpol=True for spin-polarized systems')
            elif abs(self.sigma) <= self.sigma_small and not self.fix_magmom:
                raise KeyError('Please use fix_magmom=True for sigma=0.0 eV and spinpol=True. Hopefully this is not an extended system...?')
        else:
            if self.spinpol and abs(self.sigma) <= self.sigma_small:
                self.fix_magmom = True
        if abs(self.sigma) <= self.sigma_small:
            self.occupations = 'fixed'

    def start(self):
        if not self.started:
            if self.batch:
                cdir = os.getcwd()
                os.chdir(self.localtmp)
                os.system(perHostMpiExec+' cp '+self.localtmp+'/pw.inp '+self.scratch)
                if self.calcmode!='hund':
                    self.cinp, self.cout = os.popen2(perProcMpiExec+' -wdir '+self.scratch+' pw.x -in pw.inp')
                else:
                    os.system(perProcMpiExec+' -wdir '+self.scratch+' pw.x -in pw.inp >>'+self.localtmp+'/log')
                    os.system("sed s/occupations.*/occupations=\\'fixed\\',/ <"+self.localtmp+"/pw.inp | sed s/ELECTRONS/ELECTRONS\\\\n\ \ startingwfc=\\'file\\',\\\\n\ \ startingpot=\\'file\\',/ | sed s/conv_thr.*/conv_thr="+num2str(self.conv_thr)+",/ | sed s/tot_magnetization.*/tot_magnetization="+num2str(self.totmag)+",/ >"+self.localtmp+"/pw2.inp")
                    os.system(perHostMpiExec+' cp '+self.localtmp+'/pw2.inp '+self.scratch)
                    self.cinp, self.cout = os.popen2(perProcMpiExec+' -wdir '+self.scratch+' pw.x -in pw2.inp')
                os.chdir(cdir)
            else:
                os.system('cp '+self.localtmp+'/pw.inp '+self.scratch)
                if self.calcmode!='hund':
                    self.cinp, self.cout = os.popen2('cd '+self.scratch+' ; '+'pw.x -in pw.inp')
                else:
                    os.system('cd '+self.scratch+' ; '+' pw.x -in pw.inp >>'+self.localtmp+'/log')
                    os.system("sed s/occupations.*/occupations=\\'fixed\\',/ <"+self.localtmp+"/pw.inp | sed s/ELECTRONS/ELECTRONS\\\\n\ \ startingwfc=\\'file\\',\\\\n\ \ startingpot=\\'file\\',/ | sed s/conv_thr.*/conv_thr="+num2str(self.conv_thr)+",/ | sed s/tot_magnetization.*/tot_magnetization="+num2str(self.totmag)+",/ >"+self.localtmp+"/pw2.inp")
                    os.system('cp '+self.localtmp+'/pw2.inp '+self.scratch)
                    self.cinp, self.cout = os.popen2('cd '+self.scratch+' ; '+'pw.x -in pw2.inp')

            self.started = True

    def stop(self):
        if self.started:
            if self.opt_algorithm == 'ase3':
                #sending 'Q' to espresso tells it to quit cleanly
                print >>self.cinp, 'Q'
            self.cinp.flush()
            s = open(self.localtmp+'/log','a')
            a = self.cout.readline()
            s.write(a)
            while a!='':
                a = self.cout.readline()
                s.write(a)
            s.close()
            self.cinp.close()
            self.cout.close()
            self.started = False

    def write_pot(self, filename='pot.xsf'):
        if filename[0]!='/':
            file = self.sdir+'/'+filename
        else:
            file = filename
        self.update(self.atoms)
        self.stop()
        f = open(self.localtmp+'/pp.inp', 'w')
        print >>f, '&inputPP\n  prefix=\'calc\'\n  outdir=\'.\','
        print >>f, '  plot_num=11,\n  filplot=\''+file+'\'\n/\n'
        print >>f, '&plot\n  iflag=3,\n  outputformat=3\n/'
        f.close()
        if self.batch:
            cdir = os.getcwd()
            os.chdir(self.localtmp)
            os.system(perHostMpiExec+' cp '+self.localtmp+'/pp.inp '+self.scratch)
            os.system(perProcMpiExec+' -wdir '+self.scratch+' pp.x -in pp.inp >>'+self.localtmp+'/pp.log')
            os.chdir(cdir)
        else:
            os.system('cp '+self.localtmp+'/pp.inp '+self.scratch)
            os.system('cd '+self.scratch+' ; '+'pp.x -in pp.inp >>'+self.localtmp+'/pp.log')


    def save_output(self, filename='calc.tgz'):
        if filename[0]!='/':
            file = self.sdir+'/'+filename
        else:
            file = filename
        self.update(self.atoms)
        self.stop()
        
        os.system('tar czf '+filename+' --directory='+self.scratch+' calc.save')


    def load_output(self, filename='calc.tgz'):
        self.stop()
        if filename[0]!='/':
            file = self.sdir+'/'+filename
        else:
            file = filename

        os.system('tar xzf '+filename+' --directory='+self.scratch)


    def save_wf(self, filename='wf.tgz'):
        if filename[0]!='/':
            file = self.sdir+'/'+filename
        else:
            file = filename
        self.update(self.atoms)
        self.stop()
        
        os.system('tar czf '+filename+' --directory='+self.scratch+' --exclude=calc.save .')


    def load_wf(self, filename='wf.tgz'):
        self.stop()
        if filename[0]!='/':
            file = self.sdir+'/'+filename
        else:
            file = filename

        os.system('tar xzf '+filename+' --directory='+self.scratch)


    def get_final_structure(self):
        """
        returns Atoms object according to a structure
        optimized internally by quantum espresso
        """
        from ase import Atoms

        self.stop()

        p = os.popen('grep -n Giannozzi '+self.localtmp+'/log | tail -1','r')
        n = int(p.readline().split()[0].strip(':'))
        p.close()
        
        s = open(self.localtmp+'/log','r')
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
        
        p = os.popen('grep -3 "total   stress" '+self.localtmp+'/log | tail -3','r')
        s = p.readlines()
        p.close()
        
        if len(s)!=3:
            raise RuntimeError, 'stress was not calculated\nconsider specifying calcstress or running a unit cell relaxation'
        
        stress = np.empty((3,3), np.float)
        for i in range(3):
            stress[i][:] = [float(x) for x in s[i].split()[:3]]
        
        return stress * rydberg/bohr**3


    def checkerror(self):
        p = os.popen('grep -n Giannozzi '+self.localtmp+'/log | tail -1','r')
        try:
            n = int(p.readline().split()[0].strip(':'))
        except:
            raise RuntimeError, 'Espresso executable doesn\'t seem to have been started.'
        p.close()

        p = os.popen(('tail -n +%d ' % n)+self.localtmp+'/log | grep -n %%%%%%%%%%%%%%%% |tail -2','r')
        s = p.readlines()
        p.close()

        if len(s)<2:
            return
        
        a = int(s[0].split()[0].strip(':'))+1
        b = int(s[1].split()[0].strip(':'))-a
        
        if b<1:
            return
        
        p = os.popen(('tail -n +%d ' % (a+n-1))+self.localtmp+('/log|head -%d' % b),'r')
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
            cell_factor=5.,
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
        self.oldfmax = self.fmax
        self.oldpress = self.press
        self.olddpress = self.dpress
        
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
        self.oldfmax = self.fmax
       
        self.calcmode='relax'
        if fmax is not None:
            self.fmax = fmax
        self.recalculate=True
        self.read(self.atoms)
        self.calcmode = oldmode
        self.opt_algorithm = oldalgo
        self.fmax = oldfmax
