#PBS -k oe
#PBS -l nodes=8,walltime=08:00:00
#PBS -m abe
#PBS -N sisl-test
#!/bin/sh

# mpitests/sisl-test.sh.  Generated from sisl-test.sh.in by configure.

prefix=/home/michael/Codes
exec_prefix=${prefix}

##############################################################
############### Config script for using MPICH ################
##############################################################

# get the script name
progname=$0

# Set the following entries:
JOBNAME=anyjob

# Run dir:
RUNDIR=""

# Application name:
APPLICATION="${exec_prefix}/bin/sisl-test.mpi"

# Set this for multi-processor jobs
OPFILE=""

if [[ ${OPFILE} == "" ]] ; then
  echo Unknown output file name: set OPFILE in ${progname}
  exit ;
fi

# Command line options
RUNFLAGS=""

# Extra flags for mpich:
EXTRAMPI=""

# Error trapping for unset options
if [[ ${RUNDIR} == "" ]] ; then
  echo Unknown RUNDIR: set RUNDIR in ${progname}
  exit ;
fi

if [[ ${RUNFLAGS} == "" ]] ; then
  echo Unknown RUNFLAGS: set RUNFLAGS in ${progname}
  echo You should specify an input test file name
  exit ;
fi

tmpfile=`mktemp ${JOBNAME}.XXXXXX`

RUNFLAGS="${RUNFLAGS} -o ${tmpfile}"

lockfile=.lock-sisl-test-${JOBNAME}

##############################################################
#        Below this nothing should have to be changed        #
##############################################################

echo Running from MPI $MPI_HOME
echo 
echo Changing to $RUNDIR
cd $RUNDIR

if [ -a ${lockfile} ] ; then
  echo Job $JOBNAME locked by `cat ${lockfile}`
  exit ;
fi

echo $USER > ${lockfile}

nodes=(`cat $PBS_NODEFILE`)
nnodes=${#nodes[*]}

echo Nodes: ${nodes[*]}

confile="/tmp/$PBS_JOBID.conf"

prev="" 
# Figure out whether we have 2 jobs or one determined with last entry
for i in ${nodes[*]}; do
  if [[ "$prev" != "$i" ]]; then
    pp=":1"
  else
    pp=":2"
  fi
  prev=$i
done

set prev=""
# Now create the confile
for i in ${nodes[*]} ; do
   if [[ "$i" != "$prev" ]] ; then
      echo "$i$pp" >> $confile
   fi
   prev=$i
done

echo MPICH Machinefile:
cat $confile

echo "Will run command: mpirun  -np ${nnodes} -machinefile $confile $EXTRAMPI $APPLICATION $RUNFLAGS"
echo "Starting job..."
time mpirun  -np ${nnodes} -machinefile $confile $EXTRAMPI \
    $APPLICATION $RUNFLAGS

cat ${tmpfile} > ${OPFILE}
rm ${tmpfile}

rm -f $confile
rm -f ${lockfile}
