from espresso import espresso
from sys import stderr

#keep track of ourselves so we can automatically stop us
#when a new multi-espresso object is created
espressos = []

class multiespresso:
    def __init__(self,
        ncalc = 1,
        outdirprefix = 'out',
        mtxt = 'multilog.txt',
        **kwargs
        ):
        
        #stop old espresso calculators
        while len(espressos)>0:
            espressos.pop().stop()
        
        arg = kwargs.copy()
        arg['single_calculator'] = False
        arg['numcalcs'] = ncalc
        self.ncalc = ncalc
        self.outdirprefix = outdirprefix
        self.mtxt = mtxt
        self.done = [False]*self.ncalc
        
        self.calculators = []
        for i in range(ncalc):
            arg['outdir'] = outdirprefix+'_%04d' % i
            arg['procrange'] = i
            esp = espresso(**arg)
            self.calculators.append(esp)
            espressos.append(esp)

    def wait_for_total_energies(self):
        s = open(self.mtxt, 'a')
        for i in range(self.ncalc):
            self.calculators[i].init_only(self.images[i])
            self.done[i] = False
        notdone = True
        while notdone:
            notdone = False
            for i in range(self.ncalc):
                if self.calculators[i].recalculate:
                    if not self.done[i]:
                        a = self.calculators[i].cerr.readline()
                        notdone |= (a!='' and a[:17]!='!    total energy')
                        if a[:13]=='     stopping':
                            raise RuntimeError, 'problem with calculator #%d' % i
                        elif a[:20]=='     convergence NOT':
                            raise RuntimeError, 'calculator #%d did not converge' % i
                        elif a[1:17]!='    total energy':
                            stderr.write(a)
                        else:
                            if a[0]!='!':
                                self.done[i] = False
                                print >>s, 'current free energy (calc. %3d; in scf cycle) :' % i, a.split()[-2], 'Ry'
                            else:
                                self.done[i] = True
                                print >>s, 'current free energy (calc. %3d; ionic step) :  ' % i, a.split()[-2], 'Ry'
                            s.flush()
        print >>s, ''
        s.close()
                

    def set_images(self, images):
        if len(images)!=self.ncalc:
            raise ValueError, 'number of images (%d) doesn\'t match number of calculators (%d)' % (len(images),self.ncalc)
        for i in range(self.ncalc):
            images[i].set_calculator(self.calculators[i])
        self.images = images

    def set_neb(self, neb):
        self.set_images(neb.images[1:len(neb.images)-1])
        self.neb = neb
        self.neb.neb_orig_forces = self.neb.get_forces
        self.neb.get_forces = self.nebforce

    def nebforce(self):
        self.wait_for_total_energies()
        return self.neb.neb_orig_forces()
    
    def get_world(self):
        return self.calculators[0].get_world()
