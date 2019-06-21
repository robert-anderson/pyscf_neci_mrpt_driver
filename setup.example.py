from neci_strings import write_inp
import os
from subprocess import PIPE, Popen
import itertools
import host_strings
from glob import glob
import pyscf
from pyscf.fciqmcscf import serializable_nevpt2
pyscf_dir = os.path.abspath('{}/..'.format(os.path.dirname(pyscf.__file__)))

if not os.path.exists('nevpt2_store.pkl'):  
	kwargs = {
		'atom': 'Cr 0 0 0; Cr 0 0 1.6788',
		'spin': 0,
		'basis': 'ccpvtz',
		'symmetry': True,
		'symmetry_subgroup': 'D2h'
	}
	serializable_nevpt2.SerializableNevpt2(mol_kwargs=kwargs, norb=12, nelecas=12, save_dms=True)
	assert 0

def try_mkdir(dirname):
    if os.path.exists(dirname): return not len(os.listdir(dirname))
    os.makedirs(dirname)
    return True

if not os.path.exists('FCIDUMP'):  assert 0
NELEC = int(Popen(r'grep -Eo "NELEC\s*\=\s*[0-9]+" FCIDUMP | grep -Eo "[0-9]+"',
    shell=1, stdout=PIPE, stderr=PIPE).communicate()[0])
rdm_iters = 8000
rdm_energy_iters = 100
n_nodes = 1
n_tasks_per_node = 16
reset_period = 0
time='48:00:00'
system_name = 'Cr2'
neci_exe = '/users/k1507071/code/neci_merge/build/bin/dneci'
seed = 14

def shell(cmd):
    return Popen(cmd, shell=1, stderr=PIPE, stdout=PIPE).communicate()

def link(src, dst):
    assert os.path.exists(src)
    shell('rm {1}; ln -s {0} {1}'.format(os.path.abspath(src), os.path.abspath(dst)))

def copy(src, dst):
    assert os.path.exists(src)
    shell('cp {0} {1}'.format(os.path.abspath(src), os.path.abspath(dst)))

def link_file(fname, srcdir, dstdir):
    link('{}/{}'.format(srcdir, fname), '{}/{}'.format(dstdir, fname))

def copy_file(fname, srcdir, dstdir):
    copy('{}/{}'.format(srcdir, fname), '{}/{}'.format(dstdir, fname))


def make_path(values, names):
    return os.path.abspath('/'.join(['{}_{}'.format(str(values[i]), names[i]) for i in range(len(names))]))

class ConfigIterator:
    def __init__(self, *names_and_lists):
        self.n = len(names_and_lists)//2
        self.names = [names_and_lists[2*i] for i in range(self.n)]
        self.lists = [names_and_lists[2*i+1] for i in range(self.n)]
        self.map = {self.names[i]:self.lists[i] for i in range(self.n)}
        self.it = itertools.product(*self.lists)
    def glob_path(self):
        return make_path(['*']*self.n, self.names)
    def __iter__(self):
        return self
    def next(self):
        tmp = self.it.next()
        return {self.names[i]:tmp[i] for i in range(self.n)}, make_path(tmp, self.names)

shiftdamp=0.05
memoryfacspawn=20.0
rdm_start=5000
ss_start=1500

#walkers = (int(1e5), int(5e5), int(1e6))
walkers = (int(1e6),)
coredets = (0, 20)
hwhms = (0, 0.05, 0.1, 0.25, 0.5)
seeds = (14, 15, 16, 17, 18)

confit = ConfigIterator(
        'walkers', walkers,
        'coredets', coredets,
        'hwhm', hwhms,
        'seed', seeds)

main_facs = (1,4,60,0)
spawn_facs = (2,4,50,0)
recv_facs = (2,4,50,0)
threshs=(1e-8, 1e-14)#tuple(10.0**-i for i in range(15))

for config, dirname in confit:
    if not try_mkdir(dirname): continue
    print config, dirname

    link_file('FCIDUMP', '.', dirname)
    link_file('nevpt2_store.pkl', '.', dirname)

    args = [NELEC, config['walkers']]
    kwargs = {
        'write_wf': True,
        'write_rdm': True,
        'cicj_promotion_hwhm_fraction': config['hwhm'],
        'shiftdamp': shiftdamp,
        'spawn_facs': spawn_facs, 
        'memoryfacspawn': memoryfacspawn,
        'seed': config['seed'], 
        'main_facs': main_facs, 
        'recv_facs': recv_facs,
        'rdm_iters': rdm_iters,
        'rdm_energy_iters': rdm_energy_iters
    }

    if config['coredets']==0:
        write_inp('{}/initial_hbrdm.inp'.format(dirname), *args, rdm_start=rdm_start, ss_start=ss_start, **kwargs)
        write_inp('{}/reset_hbrdm.inp'.format(dirname), *args, rdm_start=0, read_wf=True, **kwargs)
        write_inp('{}/continue_hbrdm.inp'.format(dirname), *args, read_rdm=True, rdm_start=0, **kwargs)
    else:
        write_inp('{}/initial_hbrdm.inp'.format(dirname), *args, rdm_start=rdm_start, ss_start=ss_start,
                write_core=True, pops_core=config['coredets'], **kwargs)
        write_inp('{}/reset_hbrdm.inp'.format(dirname), *args, rdm_start=0, read_wf=True, 
                read_core=True, pops_core=config['coredets'], **kwargs)
        write_inp('{}/continue_hbrdm.inp'.format(dirname), *args, read_core=True, read_rdm=True,
                rdm_start=0, **kwargs)

    kwargs = {
        'nnodes':n_nodes,
		'ntaskspernode':n_tasks_per_node,
        'name':'pt2_{}_{}'.format(system_name, '_'.join(map(str, [config[k] for k in confit.names]))),
        'time':time, 
        'wd':dirname, 
        'reset_period': reset_period,
        'threshs': str(threshs),
        'neci_exe':neci_exe,
		'pyscf_dir':pyscf_dir
    }
    #submit_script = host_strings.Nevpt2Script(pops_dir='../..')
    submit_script = host_strings.Nevpt2Script()
    submit_script.render('{}/submit.sh'.format(dirname), **kwargs)
print 'Launch command:'
print 'for i in $(ls -d {}); do cd $i; pwd; sbatch submit.sh; cd -; done'.format(confit.glob_path())
