### SGE options ################################################################
# use bash as shell
#$ -S /bin/sh
# use current (submit) directory as working directory
#$ -cwd
# join stdout and stderr
#$ -j y

# Send notification emails to:
#$ -M chrisdavis132@gmail.com
# Send on beginning, end, abortion, suspension:
#$ -m abes

# requested wall clock time
#$ -l h_rt=2:00:00

# queue
# Specify we want nodes with only 1 processor:
#$ -q single.q

# Single processor nodes are all dual core, so we need to say we require
#  8 processing units (to be allocated 4 nodes):
#$ -pe mpich 8
#$ -N davis_hpc

# Module
module load gcc/4.7.2
module load sge/6.2u5
module load openmpi/gcc/64/1.4.2

# Befehle
cd ./examples

# Run code, allocating mpi units in round-robin fashion BY NODE
#  (i.e. have one process on each node before starting a second process on a node):

# Get times for different mesh refinements:

echo "cg 5"
mpirun -np 4 --bynode main-parallel ./input/Square_four_domains 5 cg 1

echo "cg 6"
mpirun -np 4 --bynode main-parallel ./input/Square_four_domains 6 cg 1

echo "cg 7"
mpirun -np 4 --bynode main-parallel ./input/Square_four_domains 7 cg 1

echo "cg 8"
mpirun -np 4 --bynode main-parallel ./input/Square_four_domains 8 cg 1

echo "cg 9"
mpirun -np 4 --bynode main-parallel ./input/Square_four_domains 9 cg 1

echo "gs 5"
mpirun -np 4 --bynode main-parallel ./input/Square_four_domains 5 gs 1

echo "gs 6"
mpirun -np 4 --bynode main-parallel ./input/Square_four_domains 6 gs 1

echo "gs 7"
mpirun -np 4 --bynode main-parallel ./input/Square_four_domains 7 gs 1

echo "gs 8"
mpirun -np 4 --bynode main-parallel ./input/Square_four_domains 8 gs 1

echo "gs 9"
mpirun -np 4 --bynode main-parallel ./input/Square_four_domains 9 gs 1


echo "habe $NSLOTS Prozessor"
echo "Maschine:"
cat $TMPDIR/machines

exit 0