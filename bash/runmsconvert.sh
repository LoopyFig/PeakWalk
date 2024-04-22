processors=1

while getopts ":r:p:" option; do
  case $option in
    r)
      raw=$OPTARG;;
    p)
      processors=$OPTARG;;
  esac
done

cd $raw
time ( ls | xargs -I{} -P $processors perl $PEAK_WALK/perl/msconvert.pl -i {} -o mzxml -c centroid )
