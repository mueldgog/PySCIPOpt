FILE=$1

while read line
do
    COLS=( $line )
    echo ${COLS[0]}
    python run.py "data/${COLS[0]}"
done < $FILE
