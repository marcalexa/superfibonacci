#!/bin/bash

/bin/rm *ah$1.txt

../../bin/s3sample -n 19683 -o | ../../bin/s3cs -r --nvb $1 >> sfah$1.txt
../../bin/s3sample -n 19683 -o -u | ../../bin/s3cs -r --nvb $1 >> uah$1.txt
../../bin/s3cs -r --nvb $1 -f ../hopf/h19683.points >> hah$1.txt
../../bin/s3cs -r --nvb $1 -f ../soi/opt/s20360.points >> soah$1.txt
../../bin/s3cs -r --nvb $1 -f ../soi/icosa/s19188.points >> siah$1.txt
../../bin/s3cs -r --nvb $1 -f ../karney/c48u815.quat >> kah$1.txt

