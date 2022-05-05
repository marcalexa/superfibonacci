#!/bin/bash

for n in `cat ../nlist`
do
../../bin/s3sample -n $n -o | ../../bin/s3cs -r -e >> sf.txt
../../bin/s3sample -n $n -o -u | ../../bin/s3cs -r -e >> u.txt
done
for f in ../hopf/*.points
do
../../bin/s3cs -r -e -f $f >> h.txt
done
for f in ../karney/*.quat
do
../../bin/s3cs -r -e -f $f >> k.txt
done
for f in ../soi/opt/*.points
do
../../bin/s3cs -r -e -f $f >> so.txt
done
for f in ../soi/icosa/*.points
do
../../bin/s3cs -r -e -f $f >> si.txt
done

