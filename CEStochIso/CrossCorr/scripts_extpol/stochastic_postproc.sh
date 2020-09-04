#!/bin/sh -e
#source locations.sh
#source MatlabSetup_R2007a_glnxa64.sh
#./mvdata.perl .
#echo "creating links"
#./create_links.sh ${IFO}
#echo "making postDir"
#mkdir -p ${POSTDIR}
#disp=$$
#nohup /usr/X11R6/bin/Xvfb :${disp} 1>/dev/null &
#pid=$!
#echo "running post processing"
#env DISPLAY=:0.0 
echo ${MATAPPS_TOP}/src/searches/stochastic/PostProcessing/postProcessScriptFull $*
${MATAPPS_TOP}/src/searches/stochastic/PostProcessing/postProcessScriptFull $* 
#exit with same status to ensure condor stops failed jobs
exit $?
#echo "removing links"
#./rm_links.sh ${IFO}

#kill $pid
