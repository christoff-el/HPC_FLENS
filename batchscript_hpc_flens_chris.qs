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
# DUAL?????
#$ -q dual.q
#$ -pe mpich 4
#$ -N davis_hpc

# Module
#module load gcc/4.7.2
#module load sge/6.2u5
#module load openmpi/gcc/64/1.4.2

# Befehle
make
cd ./examples
mpirun -np 4 main-parallel ./input/Square_four_domains cg

echo "habe $NSLOTS Prozessor"
echo "Maschine:"
cat $TMPDIR/machines

exit 0