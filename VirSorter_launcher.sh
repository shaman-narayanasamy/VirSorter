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
VIRSORTER_REPO="/home/snarayanasamy/Work/tools/Virsorter/VirSorter"
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

CMD="docker run -v ${DATADIR}:/data -v ${VIRSORTER_REPO}/Scripts/:/usr/local/bin/Scripts \
    -v ${VIRSORTER_REPO}/wrapper_phage_contigs_sorter_iPlant.pl:/usr/local/bin/wrapper_phage_contigs_sorter_iPlant.pl \
    -v ${INPUTDIR}:/input -v ${OUTPUTDIR}:/wdir --rm discoenv/virsorter:updateBLAST \
    --db 2 --fna /input/$INPUTFILE --wdir /wdir/ --ncpu ${THREADS}"

echo $CMD
#exec $CMD
