#!/bin/bash

# check if tmp-path exists
if test ! -d /tmp/${USER}-tmpdir
then
    mkdir /tmp/${USER}-tmpdir
    echo Creating directory /tmp/${USER}-tmpdir for temporary outfile
fi

TMPFILE=/tmp/${USER}-tmpdir/$OUTFILE

python $pythonscript $INFILE $timelimit $memlimit &> $TMPFILE
cp $TMPFILE $OUTPATH/$OUTFILE
rm -f $TMPFILE
