processsors=1

while getopts ":a:p:" option; do
  case $option in
    a)
      adap=$OPTARG;;
    p)
      processors=$OPTARG;;
  esac
done

cd $adap
time ( ls | grep '.xml' | xargs -I{} -P $processors /opt/MZmine-2.53-Linux/startMZmine-Linux {} )
