#!/bin/bash


# PLEASE NOTE THIS IS AN EXAMPLE!
t1=`tconvert November 13 2016 00:00`
t3=`tconvert November 14  2016 00:00`

python make-job-file_new -s $t1 -e $t3
