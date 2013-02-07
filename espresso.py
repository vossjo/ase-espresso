#cluster-dependent definitions
scratch = '/scratch'
submitdir = '$LS_SUBCWD'
jobid = '$LSB_BATCH_JID'
getprocs = ' echo -e $LSB_HOSTS | sed s/" "/"\\\\n"/g >machinefile ;'\
          +' uniq machinefile >uniqmachinefile ;'\
          +' nodes=`wc -l <uniqmachinefile` ;'\
          +' np=`wc -l <machinefile` '
perHostMpiExec = 'mpiexec --mca plm_rsh_agent /afs/slac.stanford.edu/package/lsf/bin.slac/gmmpirun_lsgrun.sh -machinefile uniqmachinefile -np `wc -l <uniqmachinefile`'
perProcMpiExec = 'pam -g /afs/slac/g/suncat/bin/suncat-tsmpirun -x LD_LIBRARY_PATH'

from ase.calculators.general import Calculator
import atexit
import os, sys, string
import numpy as np


def checkbatch():
    p = os.popen('echo '+jobid, 'r')
    batch = (p.readline().strip()!='')
    p.close()
    return batch

def mklocaltmp(batch, odir):
    if batch:
        s = submitdir
        job = jobid
    else:
        s = '.'
        job = ''
    if odir is None:
        p = os.popen('mktemp -d '+s+'/qe"'+job+'"_XXXXX', 'r')
    else:
        p = os.popen('cd '+s+' ; mkdir -p '+odir+' ; cd '+odir+' ; pwd', 'r')
    tdir = p.readline().strip()
    p.close()
    p = os.popen('cd '+tdir+ ' ; pwd', 'r')
    tdir = p.readline().strip()
    p.close()
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
            print list2, u_l2, np.where(list2 == u_l2)
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

class espresso(Calculator):
    def __init__(self, pw=350.0, dw=3500.0, nbands=-10, 
                 kpts=(1,1,1),kptshift=(0,0,0),
                 mode='relax',
                 xc='PBE', spinpol=False,
                       outdir=None, calcstress=False,
                       psppath=None, smearing='mv', sigma=0.2,
                       U=None,J=None,U_alpha=None,U_projection_type='atomic',
                       tot_charge=0.0, # +1 means 1 e missing, -1 means 1 extra e
                       tot_magnetization=-1,#-1 means unspecified
                       occupations='smearing',#'smearing', 'fixed', 
                       dipole={'status':False},
                       field={'status':False},
                       output={'avoidio':False,
                               'removewf':True,
                               'wf_collect':False},
                       convergence={'energy':5e-6,
                                   'mixing':0.5,
                                    'maxsteps':100,
                                    'diag':'david'},
                       startingpot=None,
                       startingwfc=None):
        
        self.batch = checkbatch()
        self.localtmp = mklocaltmp(self.batch, outdir)
        if self.batch:
            self.nodes,self.np = mpisetup(self.localtmp)
        self.scratch = mkscratch(self.batch, self.localtmp)
        if output is not None and output.has_key('removewf'):
            removewf = output['removewf']
        else:
            removewf = True
        atexit.register(cleanup, self.localtmp, self.scratch, removewf, self.batch, self)
        
        #sdir is the directory the script is run or submitted from
        self.sdir = getsubmitorcurrentdir()
        
        self.pw = pw
        self.dw = dw
        self.nbands = nbands
        self.kpts = kpts
        self.kptshift = kptshift
        self.calcmode = mode
        self.xc = xc
        self.smearing = smearing
        self.sigma = sigma
        self.spinpol = spinpol
        self.tot_charge = tot_charge
        self.tot_magnetization = tot_magnetization
        self.occupations = occupations
        self.outdir = outdir
        self.calcstress = calcstress
        if psppath is None:
            try:
                self.psppath = os.environ['ESP_PSP_PATH']
            except:
                print 'Unable to find pseudopotential path.  Consider setting ESP_PSP_PATH environment variable'
                raise
        else:
            self.psppath = psppath
        if dipole is None:
            self.dipole = {'status':False}
        else:
            self.dipole = dipole
        if field is None:
            self.field = {'status':False}
        else:
            self.field = field
        self.output = output
        self.convergence = convergence
        self.startingpot = startingpot
        self.startingwfc = startingwfc
        self.U = U
        self.J = J
        self.U_alpha = U_alpha
        self.U_projection_type = U_projection_type
        self.atoms = None
        self.started = False

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
        f = open(self.localtmp+'/pw.inp', 'w')
        
        ### &CONTROL ###
        print >>f, '&CONTROL\n  calculation=\''+self.calcmode+'\',\n  prefix=\'calc\','
        print >>f, '  pseudo_dir=\''+self.psppath+'\','
        print >>f, '  outdir=\'.\','
        efield = (self.field['status']==True)
        dipfield = (self.dipole['status']==True)
        if efield or dipfield:
            print >>f, '  tefield=.true.,'
            if dipfield:
                print >>f, '  dipfield=.true.,'
        if self.U_projection_type is 'atomic':
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

        ### &SYSTEM ###
        print >>f, '/\n&SYSTEM\n  ibrav=0,\n  celldm(1)=1.8897261245650618d0,'
        print >>f, '  nat='+str(self.natoms)+','
        self.atoms2species() #self.convertmag2species()
        print >>f, '  ntyp='+str(self.nspecies)+',' #str(len(self.msym))+','
        if not self.tot_charge:
            print >>f, '  tot_charge='+str(self.tot_charge)+','
        if self.tot_magnetization != -1:
            print >>f, '  tot_magnetization='+str(self.tot_magnetization)+','

        print >>f, '  ecutwfc='+str(self.pw/rydberg)+'d0,'
        print >>f, '  ecutrho='+str(self.dw/rydberg)+'d0,'
        if self.nbands is not None:
            if self.nbands>0:
                print >>f, '  nbnd='+str(self.nbands)+','
            else:
                n = 0
                nel = {}
                for x in self.species:
                    el = x.strip('0123456789')
                    p = os.popen('grep "Z valence" '+self.psppath+'/'+el+'.UPF','r')
                    nel[el] = int(round(float(p.readline().split()[-3])))
                    p.close()
                for x in self.specprops:
                    n += nel[x[0].strip('0123456789')]
                if not self.spinpol:
                    n /= 2
                print >>f, '  nbnd='+str(n-self.nbands)+','
        if abs(self.sigma)>1e-13:
            print >>f, '  occupations=%s,'% self.occupations
            print >>f, '  smearing=\''+self.smearing+'\','
            print >>f, '  degauss='+str(self.sigma/rydberg)+'d0,'
        else:
            print >>f, '  occupations=\'fixed\','
        if self.spinpol:
            print >>f, '  nspin=2,'
            spcount  = 1
            for species in self.species: # FOLLOW SAME ORDERING ROUTINE AS FOR PSP                
                magindex = int(string.join([i for i in species if i.isdigit()],''))
                el  = species.strip('0123456789')
                mag = self.specdict[el]['magmoms'][magindex-1]
                print >>f, '  starting_magnetization(%d)=%sd0,' % (spcount,mag)
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
            print >>f, '  emaxpos='+str(emaxpos)+'d0,'
            print >>f, '  eopreg='+str(eopreg)+'d0,'
            print >>f, '  eamp='+str(eamp)+'d0,'
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
            print >>f, '  emaxpos2='+str(emaxpos)+'d0,'
            print >>f, '  eopreg2='+str(eopreg)+'d0,'
            print >>f, '  eamp2='+str(eamp)+'d0,'
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
        ### &ELECTRONS ###
        print >>f,'/\n&ELECTRONS'
        try:
            diag = self.convergence['diag']
            print >>f,'  diagonalization=\''+diag+'\','
        except:
            pass
        if self.convergence is None:
            print >>f, '  conv_thr='+str(5e-6/rydberg)+'d0,'
        else:
            if self.convergence.has_key('energy'):
                print >>f, '  conv_thr='+str(self.convergence['energy']/rydberg)+','
            else:
                print >>f, '  conv_thr='+str(5e-6/rydberg)+','
        for x in self.convergence.keys():
            if x=='mixing':
                print >>f, '  mixing_beta='+str(self.convergence[x])+'d0,'
            elif x=='maxsteps':
                print >>f, '  electron_maxstep='+str(self.convergence[x])+','
            elif x=='nmix':
                print >>f, '  mixing_ndim='+str(self.convergence[x])+','
            elif x=='mixing_mode':
                print >>f, '  mixing_mode=\''+self.convergence[x]+'\','
        if self.startingpot is not None:
            print >>f, '  startingpot=\''+self.startingpot+'\','
        if self.startingwfc is not None:
            print >>f, '  startingwfc=\''+self.startingwfc+'\','

        ### &IONS ###
        print >>f, '/\n&IONS\n  ion_dynamics=\'ase3\',\n/'

        print >>f, 'CELL_PARAMETERS'
        for i in range(3):
            print >>f, '%21.15fd0 %21.15fd0 %21.15fd0' % (self.atoms.cell[i][0],self.atoms.cell[i][1],self.atoms.cell[i][2])

        print >>f, 'ATOMIC_SPECIES'
        for species in self.species:   # PSP ORDERING FOLLOWS SPECIESINDEX
            el = species.strip('0123456789')
            print >>f, species, self.specdict[el]['mass'][0], el+'.UPF'
            print 
        
        print >>f, 'ATOMIC_POSITIONS {crystal}'
        for species, mass, magmom, pos in self.specprops:
            print >>f, '%-4s %21.15fd0 %21.15fd0 %21.15fd0' % (species,pos[0],pos[1],pos[2])
        
        print >>f, 'K_POINTS automatic'
        print >>f, self.kpts[0], self.kpts[1],self.kpts[2],self.kptshift[0],self.kptshift[1],self.kptshift[2]
        f.close()
        
    def set_atoms(self, atoms):
        if self.atoms is None:
            self.atoms = atoms.copy()
        else:
            msg = 'creation of new QE calculator object required for new atoms'
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
        if max(x.flat)>1E-13 or min(x.flat)<-1E-13 or not self.started:
            self.recalculate = True
            self.read(atoms)
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
                p = atoms.positions
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
                raise RuntimeError, 'SCF calculation failed'
            elif a=='' and self.calcmode in ('relax','scf','vc-relax','vc-md','md'):
                raise RuntimeError, 'SCF calculation didn\'t converge'
            self.atom_occ = atom_occ
            if self.calcmode in ('relax','scf','vc-relax','vc-md','md'):
                self.energy_free = float(a.split()[-2])*rydberg
            else:
                self.energy_free = None

#           a = self.cout.readline()
#           s.write(a)
#           while a[:13]!='     smearing':
#               a = self.cout.readline()
#               sys.stdout.flush()
#               s.write(a)
#           self.energy_zero = self.energy_free - float(a.split()[-2])*rydberg
            self.energy_zero = self.energy_free
            a = self.cout.readline()
            s.write(a)
            
            if self.calcmode in ('relax','scf','vc-relax','vc-md','md'):
                if self.U_projection_type == 'atomic':
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
                self.forces = None
            self.recalculate = False
            s.close()
                

    def initialize(self, atoms, calcstart=1):
        if not self.started:
            a = self.atoms
            
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
            self.writeinputfile()
        if calcstart:
            self.start()
    
    def start(self):
        if not self.started:
            if self.batch:
                cdir = os.getcwd()
                os.chdir(self.localtmp)
                self.cinp, self.cout = os.popen2(perProcMpiExec+' -wdir '+self.scratch+' pw.x -in '+self.localtmp+'/pw.inp')
                os.chdir(cdir)
            else:
                self.cinp, self.cout = os.popen2('cd '+self.scratch+' ; '+'pw.x -in '+self.localtmp+'/pw.inp')
            self.started = True

    def stop(self):
        if self.started:
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
            os.system(perProcMpiExec+' -wdir '+self.scratch+' pp.x -in '+self.localtmp+'/pp.inp >>'+self.localtmp+'/pp.log')
            os.chdir(cdir)
        else:
            os.system('cd '+self.scratch+' ; '+'pp.x -in '+self.localtmp+'/pp.inp >>'+self.localtmp+'/pp.log')


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
