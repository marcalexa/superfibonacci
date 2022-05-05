#!/bin/bash

for n in `cat ../nlist`
do
../../bin/s3sample -n $n >> sft.txt
../../bin/s3sample -n $n -u >> ut.txt
../../src/hopf/SO3_sequence $n >> ht.txt
../../bin/s3soi -n $n -f points >> sit.txt
../../bin/s3soi -r -n $n -f points >> sot.txt
done

