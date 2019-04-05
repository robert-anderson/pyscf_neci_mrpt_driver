import numpy as np
import sys, shutil, os
sys.path.append('/users/k1507071/sharedscratch/code/pyscf/')
from pyscf import gto, scf, mcscf, fciqmcscf, mrpt
import pyscf_serializer

def plot_atoms(array):
    fig = plt.figure()
    ax = Axes3D(fig)
    x_series = tuple(obj.r[0] for obj in array)
    y_series = tuple(obj.r[1] for obj in array)
    z_series = tuple(obj.r[2] for obj in array)
    ax.scatter(x_series, y_series, z_series)
    plt.show()

class Atom:
    def __init__(self, symbol, r):
        self.symbol = symbol
        self.r = np.array(r)

    def copy(self):
        return Atom(self.symbol, self.r)

    def displacement(self):
        return np.linalg.norm(self.r)

    def axis_rotate(self, theta, axis=2):
        R = np.identity(3)
        lower, upper = {0:(1,2), 1:(0,2), 2:(0,1)}[axis]
        R[lower, lower] =  np.cos(theta)
        R[lower, upper] = -np.sin(theta)
        R[upper, lower] =  np.sin(theta)
        R[upper, upper] =  np.cos(theta)
        self.r = np.dot(R, self.r)
        return self

class QuantumChemSystem(object):
    def __init__(self, verbose=4):
        self.make_atoms_list()
        self.verbose = verbose

    def make_atoms_list(self):
        self.atoms = []

    def get_mc_obj(self):
        mol = self.get_mol_obj()
        m = scf.RHF(mol)
        m.kernel()
        casscf = mcscf.CASCI(m, self.ncas, self.nelecas)
        return mol, casscf

    def make_nevpt_object(self, dirname=None):
        # dirname None implies exact solve, i.e. no NECI invokation to supply RDMs
        mol = self.get_mol_obj()
        m = scf.RHF(mol)
        m.conv_tol = 1e-10
        if dirname is not None:
            try:
                casci = mcscf.CASCI(m, self.ncas, self.nelecas)
                pyscf_serializer.selective_deserialize(casci, pfile='{}/casci.pkl'.format(dirname))
                casci.fcisolver = fciqmcscf.FCIQMCCI(mol)
                casci.fcisolver.dirname = dirname
            except IOError:
                m.kernel()
                m.analyze(with_meta_lowdin=0)

                casci = mcscf.CASCI(m, self.ncas, self.nelecas)
                casci.fcisolver = fciqmcscf.FCIQMCCI(mol)
                casci.fcisolver.only_ints = True
                try:
                    casci.kernel()
                except TypeError:
                    pass

                pyscf_serializer.selective_serialize(casci, ['mo_coeff', '_scf.mo_occ', '_scf.mo_coeff'], pfile='{}/casci.pkl'.format(dirname))
                if not os.path.abspath(dirname)==os.path.abspath('.'):
                    shutil.move('FCIDUMP', dirname)
                return
        else:
            m.kernel()
            casci = mcscf.CASCI(m, self.ncas, self.nelecas)
            casci.kernel()
        return mrpt.NEVPT(casci)

    def get_mol_obj(self):
        return None

    def serialize_cas(self, dirname='.'):
        mol, casci = self.get_mc_obj()
        casci.fcisolver = fciqmcscf.FCIQMCCI(mol)
        casci.fcisolver.only_ints = True
        try:
            casci.kernel()
        except TypeError:
            pass

        pyscf_serializer.selective_serialize(casci, ['mo_coeff', '_scf.mo_occ', '_scf.mo_coeff'], pfile='{}/casci.pkl'.format(dirname))
        if not dirname == '.':
            shutil.move('FCIDUMP', dirname)


class Homo(QuantumChemSystem):
    def __init__(self, symbol, ncas, nelecas, bond_length=1.4, basis='ccpvdz', verbose=4):
        self.symbol = symbol
        self.ncas = ncas
        self.nelecas = nelecas
        self.bond_length = bond_length
        self.basis = basis
        super(Homo, self).__init__(verbose)

    def make_atoms_list(self):
        self.atoms = []
        self.atoms.append(Atom(self.symbol, (0, 0, 0)))
        self.atoms.append(Atom(self.symbol, (0, 0, self.bond_length)))

    def get_mol_obj(self, spin_override=None):
        return gto.M(atom = [[x.symbol, tuple(x.r)] for x in self.atoms], spin=0,
                basis = {self.symbol: self.basis}, symmetry=True, symmetry_subgroup='D2h', verbose=self.verbose)

class Hetero(QuantumChemSystem):
    def __init__(self, symbols, ncas, nelecas, bond_length=1.4, verbose=4):
        self.symbols = symbols
        self.ncas = ncas
        self.nelecas = nelecas
        self.bond_length = bond_length
        super(Hetero, self).__init__(verbose)

    def make_atoms_list(self):
        self.atoms = []
        self.atoms.append(Atom(self.symbols[0], (0, 0, 0)))
        self.atoms.append(Atom(self.symbols[1], (0, 0, self.bond_length)))

    def get_mol_obj(self, spin_override=None):
        return gto.M(atom = [[x.symbol, tuple(x.r)] for x in self.atoms], spin=0,
                basis = {self.symbols[0]: 'cc-pvdz', self.symbols[1]: 'cc-pvdz'}, 
                symmetry=True, symmetry_subgroup='C2v', verbose=self.verbose, unit='angstrom')
	
class Metallocene(QuantumChemSystem):
	'''
	Metallocene geometry
    
          (cc)  Carbon separation
                _____
		       |     |
            H  C --- C  H       _
           C ----------- C       |       (mc)
         H     -- C --     H     |     Metal to aaaa
                  H              | cyclopentadienyl
                                 |
                  Fe            _|
             H         H
         H     C --- C     H
           C ----------- C     /_ (chl) angle between H bond 
               -- C --            and cyclopentadienyl plane
                  H      (ch)

	ao is the angular offset between the upper and lower Cp rings
	0 is the eclipsed conformation and pi/5 is staggered

	Fe values from 
	M. Swart / Inorganica Chimica Acta 360 (2007) 179-189:
	'''
	
	def __init__(self, metal, ncas, nelecas, mc, cc, ch, chl, spin,
			ao=0, metal_basis='ccpvdz', h_basis='ccpvdz', c_basis='ccpvdz', verbose=4):
		self.metal = metal
		self.ncas = ncas
		self.nelecas = nelecas
		self.mc = mc
		self.cc = cc
		self.ch = ch
		self.chl = chl
		self.spin = spin
		self.ao = ao
		self.metal_basis = metal_basis
		self.h_basis = h_basis
		self.c_basis = c_basis
		super(Metallocene, self).__init__(verbose)
	def make_atoms_list(self):
		self.atoms = []
		# oc is the distance from the centre of the cyclopentadienyl radical to the carbons
		# hr is the radial displacement of the hydrogen atom relative to its nearest carbon
		# hr is the axial displacement of the hydrogen atom relative to its nearest carbon
		fifth = 2*np.pi/5
		oc = 0.5*self.cc/np.sin(0.5*fifth)
		mcz = (self.mc**2-oc**2)**0.5

		hr = np.cos(np.deg2rad(self.chl))*self.ch
		hz = np.sin(np.deg2rad(self.chl))*self.ch

		atoms = []
		atoms.append(Atom(self.metal, (0, 0, 0)))

		atoms.append(Atom('C', (0, oc, mcz)))
		for i in range(4):
			atoms.append(atoms[-1].copy().axis_rotate(fifth))

		atoms.append(Atom('H', (0, oc+hr, mcz-hz)))
		for i in range(4):
			atoms.append(atoms[-1].copy().axis_rotate(fifth))

		atoms.append(Atom('C', (0, oc, -mcz)))
		atoms[-1].axis_rotate(self.ao)
		for i in range(4):
			atoms.append(atoms[-1].copy().axis_rotate(fifth))

		atoms.append(Atom('H', (0, oc+hr, -mcz+hz)))
		for i in range(4):
			atoms.append(atoms[-1].copy().axis_rotate(fifth))
		return atoms

	def get_mol_obj(self):
		return gto.M(
				atom = [[x.symbol, tuple(x.r)] for x in self.make_atoms_list()], spin=self.spin,
				basis = {self.metal: self.metal_basis, 'H': self.h_basis, 'C': self.c_basis},
				symmetry = True, verbose=self.verbose)

class Ferrocene(Metallocene):
	'''
	Fe values from 
	M. Swart / Inorganica Chimica Acta 360 (2007) 179-189:
	'''
	def __init__(self, ncas, nelecas, spin_case,
			ao=0, metal_basis='ccpvdz', h_basis='ccpvdz', c_basis='ccpvdz'):
		mc = {'low':2.007, 'interm':2.110, 'high':2.277}
		cc = {'low':1.431, 'interm':1.424, 'high':1.430}
		ch = {'low':1.087, 'interm':1.086, 'high':1.086}
		chl= {'low':1.220, 'interm':0.710, 'high':0.140}
		spin= {'low':0, 'interm':2, 'high':4}
		super(Ferrocene, self).__init__('Fe', ncas, nelecas,
				mc[spin_case], cc[spin_case], ch[spin_case], chl[spin_case],
				spin[spin_case], ao, metal_basis, h_basis, c_basis)


if __name__=='__main__':
    import time
    dirname=None
    if len(sys.argv)>1:
        dirname=os.path.abspath(sys.argv[1])
        assert os.path.exists(dirname)
    print 'UNIX time: {}'.format(time.time())
    qm_obj = Homo('N', 6, 6, 1.0977, basis='ccpvdz')
    print 'UNIX time: {}'.format(time.time())
    nevpt_obj = qm_obj.make_nevpt_object(dirname)
    print 'UNIX time: {}'.format(time.time())
    if nevpt_obj is None: sys.exit()
    nevpt_obj.kernel()
    print 'UNIX time: {}'.format(time.time())
    if dirname is not None:
        print '3RDM -> 2RDM partial trace error: {}'.format(nevpt_obj.partial_trace_error)
        print 'UNIX time: {}'.format(time.time())

