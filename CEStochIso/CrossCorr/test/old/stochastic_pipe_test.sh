#!/bin/sh
#
# This is a simple test of the stochastic pipeline code. It requires that
# the environment variable LSC_DATAFIND_SERVER be set to an appropriate
# value.
#

PARAMS_FILE=input/paramfiles/test_params.txt
JOBS_FILE=input/jobfiles/test_jobs.txt

echo "-- Removing old files"
make clean

echo
echo "-- Checking LSC_DATAFIND_SERVER"
if [ "$LSC_DATAFIND_SERVER" == "" ]; then
  echo "Error: environment variable LSC_DATAFIND_SERVER is not set."
  echo "Please set it to the appropriate value for your system."
  exit 1
fi

LSCdataFind --ping

echo
echo "-- Preparing cache stochastic pipeline"
../stochastic_pipe.tclsh -v -f -p $PARAMS_FILE -j $JOBS_FILE
res=$?
if [ $res -eq 1 ]; then
  echo "Error running stochastic_pipe.tclsh"
  exit $res
fi

echo
echo "-- Running stochastic pipeline"
if [ -x ../stochastic ]; then
  ../stochastic "$PARAMS_FILE" "$JOBS_FILE" "1"
else
  echo "../stochastic either does not exist or is not executable"
  exit 1
fi

echo
echo "-- Interactive test successful!"
echo
echo "-- To test the pipeline in Condor, use the command"
echo "     condor_submit_dag stochastic_pipe.dag"
echo

