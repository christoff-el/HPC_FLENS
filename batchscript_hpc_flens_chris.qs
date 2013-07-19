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
#$ -l h_rt=0:10:0

# queue
# Specify we want nodes with only 1 processor:
#$ -q single.q

# Single processor nodes are all dual core, so we need to say we require
#  8 processing units (to be allocated 4 nodes):
#$ -pe mpich 8
#$ -N davis_hpc

# Try to demand exclusive use of the nodes by requesting 7GB free memory
#  (nodes requested have 8GB):      (possibly unrequired)
#$ -l mem_free=7G

# Module
module load gcc/4.7.2
module load sge/6.2u5
module load openmpi/gcc/64/1.4.2

# Befehle
make
cd ./examples

# Run code, allocating mpi units in round-robin fashion BY NODE
#  (i.e. have one process on each node before starting a second process on a node):
mpirun -np 4 --bynode main-parallel ./input/Square_four_domains 7 cg 1

echo "habe $NSLOTS Prozessor"
echo "Maschine:"
cat $TMPDIR/machines

exit 0