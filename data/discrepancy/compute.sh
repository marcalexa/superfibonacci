#!/bin/bash

for n in `cat ../nlist`
do
../../bin/s3sample -n $n -o | ../../bin/s3di -r -n $1 >> sf$1.txt
../../bin/s3sample -n $n -o -u | ../../bin/s3di -r -n $1 >> u$1.txt
done
for f in ../hopf/*.points
do
../../bin/s3di -r -n $1 -f $f >> h$1.txt
done
for f in ../karney/*.quat
do
../../bin/s3di -r -n $1 -f $f >> k$1.txt
done
for f in ../soi/opt/*.points
do
../../bin/s3di -r -n $1 -f $f >> so$1.txt
done
for f in ../soi/icosa/*.points
do
../../bin/s3di -r -n $1 -f $f >> si$1.txt
done


