
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

module load general/python/2.7.10
module load libs/openblas/0.2.19/gcc5.3.0
module unload libs/openmpi/1.10.2/gcc5.3.0
module load libs/openmpi/2.0.1/intel15.0
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_core.so:$MKLROOT/lib/intel64/libmkl_sequential.so

. /users/k1507071/virtualenvs/python2/bin/activate
export OMP_NUM_THREADS=1
export PYTHONPATH=/users/k1507071/sharedscratch/code/pyscf:$PYTHONPATH

function trailing_float () {{
    echo $1 | grep -Eo "\-?[0-9]+\.[0-9]+\s*$"
}}

function get_output_value () {{
    trailing_float "$(grep -E "^$2\s" $1)"
}}

function postprocess {{
  if [ $(tail -n3 tmp.out | grep ended -c) == 0 ]; then
     echo "NECI failure detected. Exiting"
     exit
  fi
  i=$(ls neci.out* | wc -l)
  E_2rdm=$(grep "TOTAL ENERGY" tmp.out | grep -Eo "\-?[0-9]+.*")
  Herm_3rdm=$(grep "HBRDM HERMITICITY" tmp.out | grep -Eo "\-?[0-9]+.*")
  python {interface} . > pyscf.out.$i

  echo $Herm_3rdm $E_2rdm \
  $(get_output_value pyscf.out.$i Sr) \
  $(get_output_value pyscf.out.$i Si) \
  $(get_output_value pyscf.out.$i Sijrs) \
  $(get_output_value pyscf.out.$i Sijr) \
  $(get_output_value pyscf.out.$i Srsi) \
  $(get_output_value pyscf.out.$i Srs) \
  $(get_output_value pyscf.out.$i Sij) \
  $(get_output_value pyscf.out.$i Sir) \
  $(get_output_value pyscf.out.$i "Nevpt2 Energy") \
  >> results.dat

  mv tmp.out neci.out.$i
  mv FCIMCStats FCIMCStats.$i; mv FCIMCStats2 FCIMCStats2.$i
  cp *POPSFILE* last_pops
}}

if [ ! -e last_pops ]; then
  mkdir last_pops
fi

if [ $(grep "Calculation ended" neci.out.* | wc -l) = 0 ]; then
  mpirun -np {ncores} {neci_exe} initial.inp > tmp.out
  postprocess
fi

while true; do
  mpirun -np {ncores} {neci_exe} restart.inp > tmp.out
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

module unload python/3.6.3 python3/3.6 python3/recommended; module load python2/recommended
. /scratch/home/mmm0043/bin/python2/bin/activate
export OMP_NUM_THREADS=1
export PYTHONPATH=/scratch/home/mmm0043/Scratch/pyscf:$PYTHONPATH

function trailing_float () {{
    echo $1 | grep -Eo "\-?[0-9]+\.[0-9]+\s*$"
}}

function get_output_value () {{
    trailing_float "$(grep -E "^$2\s" $1)"
}}

function postprocess {{
  if [ $(tail -n3 tmp.out | grep ended -c) == 0 ]; then
     echo "NECI failure detected. Exiting"
     exit
  fi
  i=$(ls neci.out* | wc -l)
  E_2rdm=$(grep "TOTAL ENERGY" tmp.out | grep -Eo "\-?[0-9]+.*")
  Herm_3rdm=$(grep "HBRDM HERMITICITY" tmp.out | grep -Eo "\-?[0-9]+.*")
  python {interface} . > pyscf.out.$i

  echo $Herm_3rdm $E_2rdm \
  $(get_output_value pyscf.out.$i Sr) \
  $(get_output_value pyscf.out.$i Si) \
  $(get_output_value pyscf.out.$i Sijrs) \
  $(get_output_value pyscf.out.$i Sijr) \
  $(get_output_value pyscf.out.$i Srsi) \
  $(get_output_value pyscf.out.$i Srs) \
  $(get_output_value pyscf.out.$i Sij) \
  $(get_output_value pyscf.out.$i Sir) \
  $(get_output_value pyscf.out.$i "Nevpt2 Energy") \
  >> results.dat

  mv tmp.out neci.out.$i
  mv FCIMCStats FCIMCStats.$i; mv FCIMCStats2 FCIMCStats2.$i
  cp *POPSFILE* last_pops
}}

if [ ! -e last_pops ]; then
  mkdir last_pops
fi

if [ $(grep "Calculation ended" neci.out.* | wc -l) = 0 ]; then
  gerun {neci_exe} initial.inp > tmp.out
  postprocess
fi

while true; do
  gerun {neci_exe} restart.inp > tmp.out
  postprocess
done
'''}


from subprocess import Popen, PIPE
def get_jobfile():
    domainname, _ = Popen('hostname -d', shell=1, stdout=PIPE, stderr=PIPE).communicate()
    return jobfiles[domainname.strip()]
