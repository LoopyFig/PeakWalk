#!/bin/bash

processors="1"
adapdir="./adap"
featuredir="./feature"
version="4"

while getopts "hr::x::t::o::a::f::b::s::p::c::v::" option; do
  case $option in
    h)
      echo "-h Get Program Description"
      echo "-r raw directory, triggers raw processing"
      echo "-x mzxml directory, triggers mzxml processing"
      echo "  do not use -r and -x at the same time"
      echo "-t target list file, triggers target search"
      echo "-o convenient option for output directory"
      echo "  creates adap and feature subdirectories"
      echo "  do not use with -a or -f"
      echo "-a sets and creates adap directory"
      echo "  by default set to ./adap"
      echo "-f sets and creates feature directory"
      echo "  by default set to ./feature"
      echo "-b batch file for quantification, triggers quantification with -s"
      echo "  requires information on dilutions and standards"
      echo "-s specifies with standard library files to use in batch file"
      echo "-p number of parallel processors to use for processing"
      echo "-c specifies a correlation threshold for additional summary"
      echo "  requires replicate summary to be active"
      echo "-v specifies which version of mzmine to use"
      exit;;
    r)
      rawdir=$OPTARG;;
    x)
      mzxmldir=$OPTARG;;
    t)
      targetlist=$OPTARG;;
    o)
      adapdir=$OPTARG/adap/
      featuredir=$OPTARG/feature/
      mkdir -p $outputdir;;
    a)
      adapdir=$OPTARG;;
    f)
      featuredir=$OPTARG;;
    b)
      batchfile=$OPTARG;;
    s)
      stdlib=$OPTARG;;
    p)
      processors=$OPTARG;;
    c)
      correlation=$OPTARG;;
    v)
      version=$OPTARG;;
  esac
done

if [[ ! -z "$rawdir" ]]; then
  echo "converting raw"
  $PEAK_WALK/bash/runmsconvert.sh -r $rawdir -p $processors
  rawdir=$rawdir"/mzxml"
fi

if [[ ! -z "$mzxmldir" ]]; then
  rawdir=$mzxmldir
fi

if [[ ! -z "$rawdir" ]]; then
  mkdir -p $adapdir
  echo "making xml directions for extraction"
  python3 $PEAK_WALK/python/adapGC.py -r $rawdir -a $adapdir -v $version
  echo "running extraction"
  bash $PEAK_WALK/bash/runmzmine.sh -a $adapdir -p $processors -v $version
fi

if [[ ! -z "$targetlist" ]]; then
  mkdir -p $featuredir
  echo "starting targetted search"
  python3 $PEAK_WALK/python/GCCombo.py -a $adapdir -f $featuredir -t $targetlist -p $processors
  python3 $PEAK_WALK/python/GCSummary.py -f $featuredir
  if [[ ! -z "$batchfile" ]] && [[ ! -z "$stdlib" ]]; then
    python3 $PEAK_WALK/python/GCSummary.py -f $featuredir -b $batchfile -s $stdlib -x
    if [[ ! -z "$correlation" ]]; then
      python3 $PEAK_WALK/python/GCSummary.py -f $featuredir -b $batchfile -s $stdlib -x -c $correlation -n feature.subject.summary.corr.csv
    fi
  fi
fi

if [[ ! -z "$batchfile" ]] && [[ ! -z "$stdlib" ]]; then
  echo "performing quantification"
  python3 $PEAK_WALK/python/GCQuant.py -i $featuredir"/feature.sample.i.csv" -d $batchfile -s $stdlib -q $featuredir"/feature.sample.quant.csv"
  python3 $PEAK_WALK/python/GCSummary.py -f $featuredir -b $batchfile -s $stdlib -x --quant
  if [[ ! -z "$correlation" ]]; then
    python3 $PEAK_WALK/python/GCSummary.py -f $featuredir -b $batchfile -s $stdlib -x --quant -c $correlation -n feature.subject.qsummary.corr.csv
  fi
fi

