cd $1
batch=0
file=0
for f in *; do
  if [ $(( $file%40 )) == 0 ] ; then
    batch=$(( $batch+1 ))
    mkdir "batch"$batch
  fi
  file=$(( $file+1 ))
  echo $f
  mv $f "batch"$batch
done
