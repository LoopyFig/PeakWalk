#!/bin/bash

while getopts ":h" option; do
	case $option in
		h)
			echo "-h Get Program Description"
			exit;;
		r)
			$rawdir = $OPTARG;;
		t)
			$targetlist = $OPTARG;;
		o)
			$outputdir = $OPTARG
			$adapdir = $OPTARG/adap/
			$featuredir = $OPTARG/feature/
			mkdir -p $outputdir;;
		a)
			$adapdir = $OPTARG;;
		f)
			$featuredir = $OPTARG;;
		b)
			$batchfile = $OPTARG;;
		p)
			$processors = $OPTARG;;
	esac
done

mkdir -p $adapdir
mkdir -p $featuredir

python3 python/adapGC.py -i $rawdir -o $adapdir
bash bash/runmzmine.sh $adapdir $processors
python3 python/GCCombo.py $adapdir $featuredir $targetlist
python3 python/GCquant.py $featuredir/feature.best.sample.i.csv $batchfile $targetlist

