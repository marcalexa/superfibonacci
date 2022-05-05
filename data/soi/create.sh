#!/bin/bash

for n in `cat ../nlist`
do
../../bin/s3soi -n $n -f points >> sit.txt
read i < points
mv points icosa/s$i.points
../../bin/s3soi -r -n $n -f points >> sot.txt
read i < points
mv points opt/s$i.points
done

