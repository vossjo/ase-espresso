from ase.calculators.general import Calculator
from espresso import espresso
import numpy as np

class vibespresso(Calculator):
    def __init__(self,
        outdirprefix = 'out',
        **kwargs
        ):
        
        self.arg = kwargs.copy()
        self.outdirprefix = outdirprefix
        self.counter = 0
        self.equilibriumdensity = outdirprefix+'_equi.tgz'
        self.firststep = True
        self.ready = False
    
    def update(self, atoms):
        if self.atoms is not None:
            x = atoms.positions-self.atoms.positions
            if np.max(x)>1E-13 or np.min(x)<-1E-13:
                self.ready = False
        else:
            self.atoms = atoms.copy()
        self.runcalc(atoms)
        if atoms is not None:
            self.atoms = atoms.copy()
    
    def runcalc(self, atoms):
        if not self.ready:
            self.arg['outdir'] = self.outdirprefix+'_%04d' % self.counter
            self.counter += 1
            if self.firststep:
                self.esp = espresso(**self.arg)
                self.esp.set_atoms(atoms)
                self.esp.get_potential_energy(atoms)
                self.esp.save_chg(self.equilibriumdensity)
                self.firststep = False
            else:                
                self.arg['startingpot'] = 'file'
                self.esp = espresso(**self.arg)
                self.esp.set_atoms(atoms)
                self.esp.load_chg(self.equilibriumdensity)
                self.esp.get_potential_energy(atoms)
                self.esp.stop()
            self.ready = True

    def get_potential_energy(self, atoms, force_consistent=False):
        self.update(atoms)
        if force_consistent:
            return self.esp.energy_free
        else:
            return self.esp.energy_zero

    def get_forces(self, atoms):
        self.update(atoms)
        return self.esp.forces

    def get_name(self):
        return 'VibEspresso'

    def get_version(self):
        return '0.1'

