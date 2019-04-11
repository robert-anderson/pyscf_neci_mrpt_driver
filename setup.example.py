from neci_strings import write_inp
import os
from subprocess import PIPE, Popen
import itertools
import host_strings

system_def = "Homo('N', 6, 6, 1.0977, basis='ccpvdz')"
if not os.path.exists('FCIDUMP'):
    print 'python -c "import pyscf_interface as _; _.{}.make_nevpt_object().kernel(save_dms=True) > pyscf.exact.out" # do exact NEVPT2'.format(system_def)
    print 'python -c "import pyscf_interface as _; _.{}.make_nevpt_object(\'.\').kernel()" # make FCIDUMP and casci.pkl'.format(system_def)
    assert 0
NELEC = int(Popen(r'grep -Eo "NELEC\s*\=\s*[0-9]+" FCIDUMP | grep -Eo "[0-9]+"',
    shell=1, stdout=PIPE, stderr=PIPE).communicate()[0])
rdm_iters = 100
ncores_per_node = 16
time='24:00:00'
system_name = 'N2_1.0977_6o6e'
#neci_exe = '/mnt/lustre/users/k1507071/code/neci/build_gnu_E5-2650/bin/dneci'
neci_exe = '/scratch/home/mmm0043/Scratch/neci/build_gnu_release/bin/dneci'
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
        'walkers', (int(1e5),),
        'coredets', (1, 10),
        'gran', (1, 2, 4, 8),
        'seed', (14, 15, 16, 17))

main_facs = (1,1,1,0,60)
spawn_facs = (2,1,1,0,16)
recv_facs = (2,1,1,0,72)
shiftdamp=0.5
rdm_start=500
ss_start=100

for config, dirname in confit:
    print config, dirname

    os.makedirs(dirname)

    link_file('FCIDUMP', '.', dirname)
    link_file('casci.pkl', '.', dirname)

    args = [NELEC, config['walkers']]
    kwargs = {
        'write_wf': True,
        'write_rdm': True,
        'granularity': config['gran'],
        'shiftdamp': shiftdamp,
        'spawn_facs': spawn_facs, 
        'memoryfacspawn': 100.0,
        'seed': config['seed'], 
        'main_facs': main_facs, 
        'recv_facs': recv_facs, 
        'rdm_iters': rdm_iters
    }
    
    if config['coredets']==0:
        write_inp('{}/initial.inp'.format(dirname), *args, rdm_start=rdm_start, **kwargs)
        write_inp('{}/restart.inp'.format(dirname), *args, read_rdm=True, rdm_start=0, **kwargs)
    else:
        write_inp('{}/initial.inp'.format(dirname), *args, write_core=True,
                pops_core=config['coredets'], ss_start=ss_start, rdm_start=rdm_start, **kwargs)
        write_inp('{}/restart.inp'.format(dirname), *args, read_core=True, read_rdm=True,
                rdm_start=0, **kwargs)

    kwargs = {
        'ncores':1*ncores_per_node,
        'name':'pt2_{}_{}'.format(system_name, '_'.join(map(str, [config[k] for k in confit.names]))),
        'time':time, 
        'wd':dirname, 
        'neci_exe':neci_exe,
        'system_def':system_def
    }
    host_strings.render_jobfile('{}/submit.sh'.format(dirname), **kwargs)




