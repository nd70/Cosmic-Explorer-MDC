#!/bin/bash

${MATAPPS_TOP}/src/searches/stochastic/PostProcessing/combineResultsFromMultiplePairs $*
#exit with same status to ensure condor stops failed jobs
exit $?

