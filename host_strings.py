from file_composer import FileBlock
from subprocess import Popen, PIPE
import os

def get_domain():
    domainname, _ = Popen('hostname -d', shell=1, stdout=PIPE, stderr=PIPE).communicate()
    return domainname.strip()

schedulers = {'prv.rosalind.compute.estate' : 'sge', 'thomas.ucl.ac.uk' : 'sge',
    'pri.rosalind2.alces.network' : 'slurm'}
tokens = {'sge' : '#$', 'slurm' : '#SBATCH'}

class SubmitScript(FileBlock):
    def __init__(self):
        FileBlock.__init__(self)
        domain = get_domain()
        scheduler = schedulers[domain]
        token = tokens[scheduler]

        header = FileBlock()

        if scheduler=='sge':
            header <<\
            ' -S /bin/bash' <<\
            ' -N {name}' <<\
            ' -l h_rt={time}' <<\
            ' -l h_vmem=2G' <<\
            ' -M robert.anderson@kcl.ac.uk' <<\
            ' -cwd' <<\
            ' -pe mpi-16 {ncores}' <<\
            ' -m s'
        elif scheduler=='slurm':
            self << '#!/bin/bash -l'
            header <<\
            ' --nodes={nnodes}' <<\
            ' --partition={partition}' <<\
            ' --ntasks-per-node={ntaskspernode}' <<\
            ' --time={time}' <<\
            ' --job-name={name}'

        header.prefix_all(token)

        self << header
        self << 'WD={wd}'
        self << 'cd $WD'
        if scheduler=='slurm':
            self << 'if [ {nnodes} == 1 ]; then'
            self << '  tmpdir=$(mktemp -d /tmp$HOME/tmp.XXXXXXXXXX)'
            self << '  echo Working on a single node "${{SLURM_NODELIST}}", so using a temporary directory "${{tmpdir}}" on-node'
            self << '  cp FCIDUMP *.inp nevpt2_store.pkl $tmpdir'
            self << '  cd $tmpdir'
            self << 'fi'

        env_setup = FileBlock()

        if domain=='prv.rosalind.compute.estate':
            env_setup <<\
            'module purge' <<\
            'module load general/python/2.7.10' <<\
            'module load libs/openmpi/2.0.0/gcc5.3.0' <<\
            'export LD_LIBRARY_PATH=/opt/apps/libs/openblas/gcc/0.2.19/lib:$LD_LIBRARY_PATH' <<\
            'export LIBRARY_PATH=/opt/apps/libs/openblas/gcc/0.2.19/lib:$LIBRARY_PATH' <<\
            'export INCLUDE=/opt/apps/libs/openblas/gcc/0.2.19/include:$INCLUDE' <<\
            'export C_INCLUDE_PATH=/opt/apps/libs/openblas/gcc/0.2.19/include:$C_INCLUDE_PATH' <<\
            'MKLROOT=/opt/apps/intel/composer_xe_2015.3.187/mkl' <<\
            'export LD_LIBRARY_PATH=/opt/apps/intel/composer_xe_2015.3.187/mkl/lib/intel64:$LD_LIBRARY_PATH' <<\
            'export LD_LIBRARY_PATH=/opt/apps/intel/composer_xe_2015.3.187/compiler/lib/intel64:$LD_LIBRARY_PATH' <<\
            'export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_core.so:$MKLROOT/lib/intel64/libmkl_sequential.so' <<\
            '. /users/k1507071/virtualenvs/python2/bin/activate' <<\
            'export OMP_NUM_THREADS=1'
        elif domain=='thomas.ucl.ac.uk':
            env_setup <<\
            'module load general/python/2.7.10' <<\
            'module load libs/openblas/0.2.19/gcc5.3.0' <<\
            'module unload libs/openmpi/1.10.2/gcc5.3.0' <<\
            'module load libs/openmpi/2.0.1/intel15.0' <<\
            'module unload python/3.6.3 python3/3.6 python3/recommended; module load python2/recommended'
        elif domain=='pri.rosalind2.alces.network':
            env_setup <<\
            '. /users/k1507071/.bashrc' <<\
            'module purge' <<\
            'export OMP_NUM_THREADS=32' <<\
            'module load {0}/gcc/7.4.0 {0}/openmpi/3.1.4/gnu_7.4.0 {0}/openblas/0.3.6 {0}/lapack/3.8.0'.format(
                    '/users/k1507071/opt/modules')<<\
            'conda activate 2.7'

        env_setup <<\
        'export PYTHONPATH={pyscf_dir}:$PYTHONPATH' <<\
        'export PYTHONPATH={driver_dir}:$PYTHONPATH'
        self << env_setup
    
        
        processing = FileBlock()
        processing <<\
        'function copy_back () {{' <<\
        '  if [ -e $tmpdir ]; then' <<\
        '    rsync $tmpdir/*POPSFILE* $WD' <<\
        '    rsync $tmpdir/spinfree* $WD' <<\
        '    rsync $tmpdir/FCIMCStats* $WD' <<\
        '    rsync $tmpdir/CORESPACE $WD' <<\
        '    rsync $tmpdir/CICJ_HISTOGRAM $WD' <<\
        '    rsync $tmpdir/nevpt2.energy.dat $WD' <<\
        '    rsync $tmpdir/S*.dat $WD' <<\
        '    rsync $tmpdir/neci.out.* $WD' <<\
        '    rsync $tmpdir/pyscf.out.* $WD' <<\
        '    rsync $tmpdir/RDMEstimates $WD' <<\
        '  fi' <<\
        '}}' <<\
        'function clean_up_and_quit () {{' <<\
        '  copy_back' <<\
        '  exit' <<\
        '}}' <<\
        'function trailing_float () {{' <<\
        '  echo $1 | grep -Eo "\-?[0-9]+\.[0-9]+\s*$"' <<\
        '}}' <<\
        'function get_output_value () {{' <<\
        '  trailing_float "$(grep $2 $1)"' <<\
        '}}' <<\
        'function postprocess {{' <<\
        '  if [ $(tail -n3 tmp.out | grep ended -c) == 0 ]; then' <<\
        '    echo "NECI failure detected. Exiting"' <<\
        '    exit' <<\
        '  fi' <<\
        '  i=$(ls neci.out* | wc -l)' <<\
        '  E_2rdm=$(grep "TOTAL ENERGY" tmp.out | grep -Eo "\-?[0-9]+.*")' <<\
        '  Herm_3rdm=$(grep "HBRDM HERMITICITY" tmp.out | grep -Eo "\-?[0-9]+.*")' <<\
        '  read -r -d "" tmp << EOM' <<\
        'from pyscf.fciqmcscf import serializable_nevpt2'
        if os.path.exists('dms.pkl'):
            processing <<\
            'serializable_nevpt2.SerializableNevpt2(fname=\'nevpt2_store.pkl\', fciqmc_dir=\'.\', threshs={threshs}, save_dms=True)'
        else:
            processing <<\
            'serializable_nevpt2.SerializableNevpt2(fname=\'nevpt2_store.pkl\', fciqmc_dir=\'.\', threshs={threshs})'
        processing <<\
        'EOM' <<\
        '  python -c "$tmp" > pyscf.out.$i'
        
        for label in ["Sr    (-1)',", "Si    (+1)',", "Sijr  (+1) ,", "Srsi  (-1) ,",
                "Srs   (-2) ,", "Sij   (+2) ,", "Sir   (0)' ,"]:
            processing << '  grep "{0}   E" pyscf.out.$i | cut -d = -f 2 >> $(echo "{0}" | cut -d " " -f 1).energy.dat'.format(label)
            processing << '  grep "{0}   N" pyscf.out.$i | cut -d = -f 2 >> $(echo "{0}" | cut -d " " -f 1).norm.dat'.format(label)
        processing <<\
        '  grep "Nevpt2 Energy" pyscf.out.$i | cut -d = -f 2 >> nevpt2.energy.dat' <<\
        '  echo $Herm_3rdm $E_2rdm >> rdms.dat'

        if os.path.exists('dms.pkl'):
            # save stochastic RDMs at each step and compare with these deterministically obtained arrays.
            processing <<\
            '  read -r -d "" tmp << EOM' <<\
            'import numpy as np; import pickle' <<\
            'with open("dms.pkl", "rb") as f: stoch_dms=pickle.load(f)' <<\
            'with open("{}/dms.pkl", "rb") as f: exact_dms=pickle.load(f)'.format(os.getcwd())
            for rdm in (1, 2, 3, 'f3ac', 'f3ca'):
                processing << 'print 1.0-np.linalg.norm(stoch_dms[\'{0}\'])/np.linalg.norm(exact_dms[\'{0}\']),'.format(rdm)
            processing <<\
            'EOM' <<\
            '  python -c "$tmp" >> rdm_errors.dat'

        processing <<\
        '  if [ -e Si.pkl ]; then' <<\
        '    for contraction in $(ls S*.pkl); do' <<\
        '    mv $contraction contractions/$contraction.$i' <<\
        '    done' <<\
        '  fi' <<\
        '  mv tmp.out neci.out.$i' <<\
        '  mv FCIMCStats FCIMCStats.$i; mv FCIMCStats2 FCIMCStats2.$i' <<\
        '  cp *POPSFILE* last_pops' <<\
        '}}'

        processing <<\
        'if [ ! -e last_pops ]; then' <<\
        '  mkdir last_pops' <<\
        'fi' <<\
        'if [ ! -e contractions ]; then' <<\
        '  mkdir contractions' <<\
        'fi'
        self << processing
    def render(self, fname, **kwargs):
        if 'partition' not in kwargs: kwargs['partition'] = '"morty,nodes"'
        kwargs['driver_dir'] = os.path.abspath(__file__)
        kwargs['ncores'] = kwargs['nnodes']*kwargs['ntaskspernode']
        FileBlock.render(self, fname, **kwargs)

class WfInitScript(SubmitScript):
    def __init__(self):
        SubmitScript.__init__(self)
        execution = FileBlock()
        execution <<\
        'echo "initialising FCIQMC calculation"' <<\
        'mpirun -np {ncores} {neci_exe} low_rdm.inp > tmp.out' <<\
        'clean_up_and_quit'
        self << execution

class Nevpt2Script(SubmitScript):
    def __init__(self, pops_dir=None, step_seed=False, loop=False):
        SubmitScript.__init__(self)
        if pops_dir is not None: self << 'cp {}/*POPSFILE* .'.format(pops_dir)
        if pops_dir is not None: self << 'cp {}/CORESPACE .'.format(pops_dir)
        execution = FileBlock()
        execution <<\
        'while true; do' <<\
        '  if [ -e $WD/stop ]; then' <<\
        '    rm $WD/stop' <<\
        '    clean_up_and_quit' <<\
        '  fi' <<\
        '  i=$(ls neci.out* | wc -l)' <<\
        '  if [ $i = 0 ]; then' <<\
        '    echo "initialising FCIQMC calculation"' <<\
        '    mpirun -np {ncores} {neci_exe} initial_hbrdm.inp > tmp.out' <<\
        '  elif [ {reset_period} != 0 ] && [ $(echo "$i%{reset_period}" | bc) = 0 ]; then' <<\
        '    echo "resetting RDM accumulations and continuing from scratch"' <<\
        '    mpirun -np {ncores} {neci_exe} reset_hbrdm.inp > tmp.out' <<\
        '  else'
        if step_seed: execution << '    sed -i s/"seed.*"/"seed $i"/g continue_hbrdm.inp'
        execution <<\
        '    echo "continuing an existing RDM accumulation"' <<\
        '    mpirun -np {ncores} {neci_exe} continue_hbrdm.inp > tmp.out' <<\
        '  fi' <<\
        '  postprocess'
        if not loop: execution << '  clean_up_and_quit'
        execution << '  save_pops'
        execution << 'done'
        self << execution
