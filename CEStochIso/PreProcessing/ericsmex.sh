#! /bin/bash
# This compilation is adapted from the old compilation, which used to occur
# in extras/ using ET_mkframeYoussef2.c, now renamed to be mkframe_preproc.c.
DIR=/archive/home/ethrane/FrContrib/extras
mex mkframe_preproc.c $DIR/../v6r24/src/libFrame.a -I$DIR/../v6r24/src
mv mkframe_preproc.mexa64 ./mkframe.mexa64

# To compile the necessary libraries
#FRVER=v6r24
#gcc -O -fexceptions -fPIC extras/stat.c ${FRVER}/src/libFrame.a -I${FRVER}/src -lm -o extras/stat
