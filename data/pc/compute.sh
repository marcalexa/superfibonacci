#!/bin/bash


../../bin/s3sample -n 19683 -o | ../../bin/s3pc -a -r -n 10000 -b 1000 >> sfpc.txt
../../bin/s3sample -n 19683 -o -u | ../../bin/s3pc -a -r -n 10000 -b 1000 >> upc.txt
../../bin/s3pc -a -r -n 10000 -b 1000 -f ../hopf/h19683.points >> hpc.txt
../../bin/s3pc -a -r -n 10000 -b 1000 -f ../soi/opt/s20360.points >> sopc.txt
../../bin/s3pc -a -r -n 10000 -b 1000 -f ../soi/icosa/s19188.points >> sipc.txt
../../bin/s3pc -a -r -n 10000 -b 1000 -f ../karney/c48u815.quat >> kpc.txt

