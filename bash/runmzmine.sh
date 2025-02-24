processsors=1
mzmine2="/opt/MZmine-2.53-Linux/startMZmine-Linux {}"
mzmine4="mzmine4 -batch {}"
mzmine="$mzmine2"

while getopts ":a:p:v:" option; do
  case $option in
    a)
      adap=$OPTARG;;
    p)
      processors=$OPTARG;;
    v)
      case $OPTARG in
        2)
          mzmine=$mzmine2;;
        4)
          mzmine=$mzmine4;;
      esac
  esac
done

cd $adap
time ( ls | grep '.xml' | xargs -I{} -P $processors $mzmine )
