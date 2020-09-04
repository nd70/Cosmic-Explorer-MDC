#! /bin/bash

#This version of mkframe was modified by Bernard's student Youssef.
#This version is necessary to write complex data to SID frames.
#To compile it, check out: lscsoft/ligotools/packages/FrContrib
#Run this script to compile mkframe_SID.c
#Make sure that matlab points to it before you run/build preproc.
#(To check which version you are using, type "which mkframe".)
#
# Eric Thrane

FrDir=~/lscsoft/ligotools/packages/FrContrib
FRVER=v6r24
#

if [ ! -e $FrDir/extras/mkframe_SID.c ]
then
    echo "Can not find mkframe_SID.c in $FrDir."
    echo "Please check out FrContrib."
    exit
fi

if [ ! -e $FrDir/$FRVER ]
then
    read -p "Need to build library, ^C to cancel..."
    here=`pwd`
    cd $FrDir
    gcc -O -fexceptions -fPIC extras/stat.c ${FRVER}/src/libFrame.a -I${FRVER}/src -lm -o extras/stat   
    cd $here
fi

mex $FrDir/extras/mkframe_SID.c $FrDir/$FRVER/src/libFrame.a -I$FrDir/$FRVER/src
mex $FrDir/extras/mkframe_SID.c $FrDir/$FRVER/src/libFrame.a -I$FrDir/$FRVER/src
mv mkframe_SID.mexglx mkframe.mexglx

mex $FrDir/extras/frgetvect_SID.c $FrDir/$FRVER/src/libFrame.a -I$FrDir/$FRVER/src
mex $FrDir/extras/frgetvect_SID.c $FrDir/$FRVER/src/libFrame.a -I$FrDir/$FRVER/src
mv frgetvect_SID.mexglx frgetvect.mexglx

mex $FrDir/extras/frgetparams_SID.c $FrDir/$FRVER/src/libFrame.a -I$FrDir/$FRVER/src
mex $FrDir/extras/frgetparams_SID.c $FrDir/$FRVER/src/libFrame.a -I$FrDir/$FRVER/src
mv frgetparams_SID.mexglx frgetparams.mexglx

