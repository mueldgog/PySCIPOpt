#!/bin/bash

if [[ $# -lt 7 ]]; then
   echo "Need arguments to run this script."
   echo "usage: ./runjobs.sh <python script> <test set> <output path> <time limit> <memory limit> <queue> <exclusive: T or F>"
   exit 0
fi

# collect input data
pythonscript=$1
testset=$2
outputpath=$3
timelimit=$4
memlimit=$5
queue=$6

# default parameters
currentpath=`pwd`
exclusive=""

# reinterpret exclusive flag
if [[ $7 == "T" ]]; then
    exclusive=" --exclusive"
fi

# submit a job for each instance of the test set
while read file
do
    # extract instance name
    instancename=(${file//\// })

    INFILE=$currentpath"/data/"$instancename
    OUTPATH=$currentpath"/"$outputpath
    OUTFILE=$USER.$pythonscript.$instancename".out"

    export pythonscript
    export INFILE
    export OUTPATH
    export OUTFILE
    export timelimit
    export memlimit

    # send the job to the cluster
    if [[ $queue == "mip-dbg" ]]; then
        sbatch -A "mip" -p ${queue} --time=${timelimit} --output=/dev/null ./runjob.sh
    else
        sbatch -A "mip" -p ${queue} --time=${timelimit} ${exclusive} --mem=${memlimit} --output=/dev/null ./runjob.sh
    fi
done < $testset
