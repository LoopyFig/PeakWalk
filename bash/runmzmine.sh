cd $1
time ( ls | grep '.xml' | xargs -I{} -P$2 /opt/MZmine-2.53-Linux/startMZmine-Linux {} )
