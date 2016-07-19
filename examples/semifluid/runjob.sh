#!/bin/bash

# echo $pythonscript
# echo $INFILE
# echo $OUTFILE
# echo $timelimit
# echo $memlimit

# echo "python "$pythonscript $INFILE $timelimit $memlimit " &> " $OUTFILE

python $pythonscript $INFILE $timelimit $memlimit &> $OUTFILE
