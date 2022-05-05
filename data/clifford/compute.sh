#!/bin/bash

#../../bin/s3ct -r -t 0.1 -z -f ../hopf/h1048576.points -s 0 > h1M.pgm
#../../bin/s3ct -r -t 0.1 -z -f ../soi/opt/s1052834.points -s 0 > so1M.pgm
#../../bin/s3ct -r -t 0.1 -z -f ../soi/icosa/s1115676.points -s 0 > si1M.pgm
#../../bin/s3sample -n 1050000 -o | ../../bin/s3ct -r -t 0.1 -z -s 0 > sf1M.pgm
#../../bin/s3sample -n 1050000 -o -u | ../../bin/s3ct -r -t 0.1 -z -s 0 > u1M.pgm
#
#../../bin/s3ct -r -t 0.1 -z -f ../hopf/h514229.points -s 0 > h500K.pgm
#../../bin/s3ct -r -t 0.1 -z -f ../soi/opt/s522976.points -s 0 > so500K.pgm
#../../bin/s3ct -r -t 0.1 -z -f ../soi/icosa/s525028.points -s 0 > si500K.pgm
#../../bin/s3sample -n 525000 -o | ../../bin/s3ct -r -t 0.1 -z -s 0 > sf500K.pgm
#../../bin/s3sample -n 525000 -o -u | ../../bin/s3ct -r -t 0.1 -z -s 0 > u500K.pgm
#
#../../bin/s3ct -r -t 0.1 -z -f ../hopf/h1048576.points -s 1 > h1Mr.pgm
#../../bin/s3ct -r -t 0.1 -z -f ../soi/opt/s1052834.points -s 1 > so1Mr.pgm
#../../bin/s3ct -r -t 0.1 -z -f ../soi/icosa/s1115676.points -s 1 > si1Mr.pgm
#../../bin/s3sample -n 1050000 -o | ../../bin/s3ct -r -t 0.1 -z -s 1 > sf1Mr.pgm
#../../bin/s3sample -n 1050000 -o -u | ../../bin/s3ct -r -t 0.1 -z -s 1 > u1Mr.pgm

#../../bin/s3ct -r -t 0.1 -z -f ../hopf/h514229.points -s 1 > h500Kr.pgm
#../../bin/s3ct -r -t 0.1 -z -f ../soi/opt/s522976.points -s 1 > so500Kr.pgm
#../../bin/s3ct -r -t 0.1 -z -f ../soi/icosa/s525028.points -s 1 > si500Kr.pgm
#../../bin/s3sample -n 525000 -o | ../../bin/s3ct -r -t 0.1 -z -s 1 > sf500Kr.pgm
#../../bin/s3sample -n 525000 -o -u | ../../bin/s3ct -r -t 0.1 -z -s 1 > u500Kr.pgm
#
#../../bin/s3sample -n 2000000 -o | ../../bin/s3ct -r -t 0.1 -z -s 0 > sf2M.pgm
#../../bin/s3sample -n 2000000 -o -u | ../../bin/s3ct -r -t 0.1 -z -s 0 > u2M.pgm
#
#../../bin/s3sample -n 4000000 -o | ../../bin/s3ct -r -t 0.1 -z -s 0 > sf4M.pgm
#../../bin/s3sample -n 4000000 -o -u | ../../bin/s3ct -r -t 0.1 -z -s 0 > u4M.pgm
#
#../../bin/s3sample -n 8000000 -o | ../../bin/s3ct -r -t 0.1 -z -s 0 > sf8M.pgm
#../../bin/s3sample -n 8000000 -o -u | ../../bin/s3ct -r -t 0.1 -z -s 0 > u8M.pgm
#
#../../bin/s3sample -n 16000000 -o | ../../bin/s3ct -r -t 0.1 -z -s 0 > sf16M.pgm
#../../bin/s3sample -n 16000000 -o -u | ../../bin/s3ct -r -t 0.1 -z -s 0 > u16M.pgm

#../../bin/s3ct -r -t 0.1 -z -f ../hopf/h32768.points -s 0 > h25msec.pgm
#../../bin/s3ct -r -t 0.1 -z -f ../soi/opt/s1065.points -s 0 > so25msec.pgm
#../../bin/s3ct -r -t 0.1 -z -f ../soi/icosa/s525028.points -s 0 > si25msec.pgm
../../bin/s3sample -n 725000 -o | ../../bin/s3ct -r -t 0.1 -z -s 0 > sf25msec.pgm
../../bin/s3sample -n 300000 -o -u | ../../bin/s3ct -r -t 0.1 -z -s 0 > u25msec.pgm



for f in *.pgm
do
sips -s format png $f --out ${f%.pgm}.png
rm $f
done

