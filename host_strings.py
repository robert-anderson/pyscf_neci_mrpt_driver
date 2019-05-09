
jobfiles={
    'prv.rosalind.compute.estate':\
r'''#$ -S /bin/bash
#$ -N {name}
#$ -l h_rt={time}
#$ -l h_vmem=2G
#$ -M robert.anderson@kcl.ac.uk
#$ -cwd
#$ -pe mpi-16 {ncores}
#$ -m s

cd {wd}

{env_setup}
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_core.so:$MKLROOT/lib/intel64/libmkl_sequential.so

. /users/k1507071/virtualenvs/python2/bin/activate
export OMP_NUM_THREADS=1
export PYTHONPATH=/users/k1507071/sharedscratch/code/pyscf:$PYTHONPATH
export PYTHONPATH={driver_dir}:$PYTHONPATH

function trailing_float () {{
    echo $1 | grep -Eo "\-?[0-9]+\.[0-9]+\s*$"
}}

function get_output_value () {{
    trailing_float "$(grep $2 $1)"
}}

function postprocess {{
  if [ $(tail -n3 tmp.out | grep ended -c) == 0 ]; then
     echo "NECI failure detected. Exiting"
     exit
  fi
  i=$(ls neci.out* | wc -l)
  E_2rdm=$(grep "TOTAL ENERGY" tmp.out | grep -Eo "\-?[0-9]+.*")
  Herm_3rdm=$(grep "HBRDM HERMITICITY" tmp.out | grep -Eo "\-?[0-9]+.*")
  python -c "from pyscf.fciqmcscf import serializable_nevpt2 as _; _.SerializableNevpt2(fname='nevpt2_store.pkl', fciqmc_dir='.', threshs={threshs})" > pyscf.out.$i

  {result_output}
  {rdm_errors}

  for contraction in $(ls S*.pkl); do
    mv $contraction contractions/$contraction.$i
  done
  mv tmp.out neci.out.$i
  mv FCIMCStats FCIMCStats.$i; mv FCIMCStats2 FCIMCStats2.$i
  cp *POPSFILE* last_pops
}}

if [ ! -e last_pops ]; then
  mkdir last_pops
fi

if [ ! -e contractions ]; then
  mkdir contractions
fi

while true; do
  i=$(ls neci.out* | wc -l)
  if [ $i = 0 ]; then
    echo "initialising FCIQMC calculation"
    mpirun -np {ncores} {neci_exe} initial_hbrdm.inp > tmp.out
  elif [ {reset_period} != 0 ] && [ $(echo "$i%{reset_period}" | bc) = 0 ]; then
    echo "resetting RDM accumulations and continuing from scratch"
    mpirun -np {ncores} {neci_exe} reset_hbrdm.inp > tmp.out
  else
    echo "continuing an existing RDM accumulation"
    mpirun -np {ncores} {neci_exe} continue_hbrdm.inp > tmp.out
  fi
  postprocess
done
''','thomas.ucl.ac.uk':\
r'''#$ -S /bin/bash
#$ -N {name}
#$ -l h_rt={time}
#$ -l mem=4G
#$ -M robert.anderson@kcl.ac.uk
#$ -cwd
#$ -pe mpi {ncores}
#$ -P Gold
#$ -A KCL_Booth
#$ -m s

cd {wd}

{env_setup}
. /scratch/home/mmm0043/bin/python2/bin/activate
export OMP_NUM_THREADS=1
export PYTHONPATH=/scratch/home/mmm0043/Scratch/pyscf:$PYTHONPATH
export PYTHONPATH={driver_dir}:$PYTHONPATH

function trailing_float () {{
    echo $1 | grep -Eo "\-?[0-9]+\.[0-9]+\s*$"
}}

function get_output_value () {{
    trailing_float "$(grep $2 $1)"
}}

function postprocess {{
  if [ $(tail -n3 tmp.out | grep ended -c) == 0 ]; then
     echo "NECI failure detected. Exiting"
     exit
  fi
  i=$(ls neci.out* | wc -l)
  E_2rdm=$(grep "TOTAL ENERGY" tmp.out | grep -Eo "\-?[0-9]+.*")
  Herm_3rdm=$(grep "HBRDM HERMITICITY" tmp.out | grep -Eo "\-?[0-9]+.*")
  python -c "from pyscf.fciqmcscf import serializable_nevpt2 as _; _.SerializableNevpt2(fname='nevpt2_store.pkl', fciqmc_dir='.', threshs={threshs})" > pyscf.out.$i

  {result_output}
  {rdm_errors}

  for contraction in $(ls S*.pkl); do
    mv $contraction contractions/$contraction.$i
  done
  mv tmp.out neci.out.$i
  mv FCIMCStats FCIMCStats.$i; mv FCIMCStats2 FCIMCStats2.$i
  cp *POPSFILE* last_pops
}}

if [ ! -e last_pops ]; then
  mkdir last_pops
fi

if [ ! -e contractions ]; then
  mkdir contractions
fi

while true; do
  i=$(ls neci.out* | wc -l)
  if [ $i = 0 ]; then
    echo "initialising FCIQMC calculation"
    gerun {neci_exe} initial_hbrdm.inp > tmp.out
  elif [ {reset_period} != 0 ] && [ $(echo "$i%{reset_period}" | bc) = 0 ]; then
    echo "resetting RDM accumulations and continuing from scratch"
    gerun {neci_exe} reset_hbrdm.inp > tmp.out
  else
    echo "continuing an existing RDM accumulation"
    gerun {neci_exe} continue_hbrdm.inp > tmp.out
  fi
  postprocess
done
'''}



#for label in "Sr    (-1)'," "Si    (+1)'," "Sijrs (0)  ," "Sijr  (+1) ," "Srsi  (-1) ," "Srs   (-2) ," "Sij   (+2) ," "Sir   (0)' ,"; do
# common blocks:
result_output = \
r'''
  for label in "Sr    (-1)'," "Si    (+1)'," "Sijr  (+1) ," "Srsi  (-1) ," "Srs   (-2) ," "Sij   (+2) ," "Sir   (0)' ,"; do
    grep "$label   E" pyscf.out.$i | cut -d = -f 2 >> $(echo $label | cut -d " " -f 1).energy.dat
    grep "$label   N" pyscf.out.$i | cut -d = -f 2 >> $(echo $label | cut -d " " -f 1).norm.dat
  done
  grep "Nevpt2 Energy" pyscf.out.$i | cut -d = -f 2 >> nevpt2.energy.dat
  echo $Herm_3rdm $E_2rdm >> rdms.dat
'''

from subprocess import Popen, PIPE
import os

def output_rdm_errors():
    # assumes the exact dms are in the current directory
    python_commands = []
    python_commands.append('import numpy as np')
    python_commands.append('import pickle')
    python_commands.append('f = open(\'dms.pkl\', \'rb\')')
    python_commands.append('stoch_dms=pickle.load(f)')
    python_commands.append('f.close()')
    python_commands.append('f = open(\'{}/dms.pkl\', \'rb\')'.format(os.getcwd()))
    python_commands.append('exact_dms=pickle.load(f)')
    python_commands.append('f.close()')
    python_commands.append('print \' \'.join([str(1.0-np.linalg.norm(stoch_dms[i])/np.linalg.norm(exact_dms[i])) for i in (\'1\', \'2\', \'3\', \'f3ac\', \'f3ca\')])')
    return r'python -c "{}" >> rdm_errors.dat'.format('; '.join(python_commands))

def get_domain():
    domainname, _ = Popen('hostname -d', shell=1, stdout=PIPE, stderr=PIPE).communicate()
    return domainname.strip()

def get_jobfile():
    return jobfiles[get_domain()]

def render_jobfile(fname, **kwargs):
    kwargs['driver_dir'] = os.path.abspath(os.path.dirname(__file__))
    kwargs['result_output'] = result_output
    # if the saved exact RDMs are present, assume that we want the errors calculated at each postprocessing
    if os.path.exists('dms.pkl'):
        kwargs['save_dms'] = True
        kwargs['rdm_errors'] = output_rdm_errors()
    else:
        kwargs['save_dms'] = False
        kwargs['rdm_errors'] = ''

    if get_domain()=='prv.rosalind.compute.estate':
        if 'gnu' in kwargs['neci_exe']:
            kwargs['env_setup'] = \
'''
module purge
module load general/python/2.7.10 
module load libs/openmpi/2.0.0/gcc5.3.0
export LD_LIBRARY_PATH=/opt/apps/libs/openblas/gcc/0.2.19/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=/opt/apps/libs/openblas/gcc/0.2.19/lib:$LIBRARY_PATH
export INCLUDE=/opt/apps/libs/openblas/gcc/0.2.19/include:$INCLUDE
export C_INCLUDE_PATH=/opt/apps/libs/openblas/gcc/0.2.19/include:$C_INCLUDE_PATH
MKLROOT=/opt/apps/intel/composer_xe_2015.3.187/mkl
export LD_LIBRARY_PATH=/opt/apps/intel/composer_xe_2015.3.187/mkl/lib/intel64:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/opt/apps/intel/composer_xe_2015.3.187/compiler/lib/intel64:$LD_LIBRARY_PATH
'''
        else:
            kwargs['env_setup'] = \
'''
module load general/python/2.7.10
module load libs/openblas/0.2.19/gcc5.3.0
module unload libs/openmpi/1.10.2/gcc5.3.0
module load libs/openmpi/2.0.1/intel15.0
'''
    else:
        kwargs['env_setup'] = \
'''            
module unload python/3.6.3 python3/3.6 python3/recommended; module load python2/recommended
'''
    with open(fname, 'w') as f: f.write(get_jobfile().format(**kwargs))

