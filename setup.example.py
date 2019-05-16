from neci_strings import write_inp
import os
from subprocess import PIPE, Popen
import itertools
import host_strings
from glob import glob

if not os.path.exists('FCIDUMP'):  assert 0
NELEC = int(Popen(r'grep -Eo "NELEC\s*\=\s*[0-9]+" FCIDUMP | grep -Eo "[0-9]+"',
    shell=1, stdout=PIPE, stderr=PIPE).communicate()[0])
rdm_iters = 100
ncores_per_node = 16
reset_period = 20
time='24:00:00'
system_name = 'Cr2_1.6788_12o12e'
neci_exe = '/mnt/lustre/users/k1507071/code/neci/build_gnu_E5-2650/bin/dneci'
#neci_exe = '/scratch/home/mmm0043/Scratch/neci/build_gnu_release/bin/dneci'
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


class ConfigIterator:
    def __init__(self, *names_and_lists):
        self.n = len(names_and_lists)//2
        self.names = [names_and_lists[2*i] for i in range(self.n)]
        self.lists = [names_and_lists[2*i+1] for i in range(self.n)]
        self.map = {self.names[i]:self.lists[i] for i in range(self.n)}
        self.it = itertools.product(*self.lists)
    def __iter__(self):
        return self
    def next(self):
        tmp = self.it.next()
        return {self.names[i]:tmp[i] for i in range(self.n)}, os.path.abspath(
                '/'.join(['{}_{}'.format(str(tmp[i]), self.names[i]) for i in range(self.n)]))

confit = ConfigIterator(
        'walkers', (int(1e6), ),
        'coredets', (20, 200),
        'gran', (8, 16, 32),
        'seed', (14, 15, 16))


main_facs = (1,1,4,0,60)
spawn_facs = (2,1,4,0,50)
recv_facs = (2,1,8,0,50)
shiftdamp=0.05
rdm_start=1500
ss_start=1000
memoryfacspawn=20.0
threshs=tuple(10.0**-i for i in range(15))

for config, dirname in confit:
    print config, dirname

    os.makedirs(dirname)

    link_file('FCIDUMP', '.', dirname)
    link_file('nevpt2_store.pkl', '.', dirname)

    args = [NELEC, config['walkers']]
    kwargs = {
        'write_wf': True,
        'write_rdm': True,
        'granularity': config['gran'],
        'shiftdamp': shiftdamp,
        'spawn_facs': spawn_facs, 
        'memoryfacspawn': memoryfacspawn,
        'seed': config['seed'], 
        'main_facs': main_facs, 
        'recv_facs': recv_facs,
        'rdm_iters': rdm_iters
    }

    
    if config['coredets']==0:
        write_inp('{}/low_rdm.inp'.format(dirname), *args, two_rdm_only=True, rdm_start=rdm_start, **kwargs)
        write_inp('{}/initial_hbrdm.inp'.format(dirname), *args, rdm_start=rdm_start, **kwargs)
        write_inp('{}/reset_hbrdm.inp'.format(dirname), *args, rdm_start=0, read_wf=True, **kwargs)
        write_inp('{}/continue_hbrdm.inp'.format(dirname), *args, read_rdm=True, rdm_start=0, **kwargs)
    else:
        write_inp('{}/low_rdm.inp'.format(dirname), *args, two_rdm_only=True, rdm_start=rdm_start, 
                write_core=True, pops_core=config['coredets'], ss_start=ss_start, **kwargs)
        write_inp('{}/initial_hbrdm.inp'.format(dirname), *args, rdm_start=rdm_start, 
                write_core=True, pops_core=config['coredets'], ss_start=ss_start, **kwargs)
        write_inp('{}/reset_hbrdm.inp'.format(dirname), *args, rdm_start=0, read_wf=True, 
                read_core=True, pops_core=config['coredets'], **kwargs)
        write_inp('{}/continue_hbrdm.inp'.format(dirname), *args, read_core=True, read_rdm=True,
                rdm_start=0, **kwargs)

    kwargs = {
        'ncores':1*ncores_per_node,
        'name':'pt2_{}_{}'.format(system_name, '_'.join(map(str, [config[k] for k in confit.names]))),
        'time':time, 
        'wd':dirname, 
        'reset_period': reset_period,
        'threshs': str(threshs),
        'neci_exe':neci_exe
    }
	submit_script = host_strings.SubmitScript()
    submit_script.render('{}/submit.sh'.format(dirname), **kwargs)




