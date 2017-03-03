#!/bin/bash -l


display_usage() { 
	echo "This script is a wrapper for launching VirSorter"
	echo "Usage: $0 <input file> <output directory> <threads>"
	}

	if [ $# -le 2 ]
	then
	    display_usage
	    exit 1
	fi

DATADIR="/home/snarayanasamy/Work/tools/Virsorter/virsorter-data"
FILE=$1
OUTPUTDIR=$2
THREADS=$3

INPUTFILE=`echo ${FILE##*/}`
INPUTDIR=`echo $FILE | sed -e "s:/${INPUTFILE}::g"`

echo $INPUTFILE
echo $INPUTDIR

if [ -d "$OUTPUTDIR" ]; then
    echo "The output directory; $OUTPUTDIR already exists"
    echo "Proceeding to launch VirSorter" 
else
    echo "The output directory; $OUTPUTDIR does not exist"
    echo "Creating the output directory"
    mkdir $OUTPUTDIR
fi

#$ docker run -v /home/snarayanasamy/Work/tools/Virsorter/virsorter-data:/data
#-v /home/snarayanasamy/Work/virsorter_test/input:/input -v
#/home/snarayanasamy/Work/virsorter_test/output:/root -v
#/home/snarayanasamy/Work/virsorter_test/intermediate:/wdir --rm
#discoenv/virsorter:v1.0.3 --db 2 --fna /input/Bio17-1.fa --wdir /wdir/

CMD="docker run -v ${DATADIR}:/data -v ${INPUTDIR}:/indir -v ${OUTPUTDIR}:/wdir \
    --rm discoenv/virsorter:v1.0.3 -w /wdir \
    --db 2 --fna /indir/${INPUTFILE} \
    --cpu ${THREADS}"
echo $CMD
#exec $CMD
