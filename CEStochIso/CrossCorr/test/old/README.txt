This is the test directory for the stochastic pipeline. Before running,
please make sure that

1) the environment variable LSC_DATAFIND_SERVER is set to the appropriate
value for your local system

2) stochastic.m has been compiled

(Alternatively, you could modify the test script so that it just produces
the cache files, then run stochastic.m interactively in Matlab with the
appropriate parameters)

To run the test, run the shell script
   
   ./stochastic_pipe_test.sh

This will run stochastic_pipe.tclsh to generate the cache files in
input/cachefiles, then run the stochastic executable, putting the results
into the output/ directory.

To test the pipeline under condor, use the command

   condor_submit_dag stochastic_pipe.dag

after running stochastic_pipe_test.sh.

The clean.sh script can be used to remove the output files produced by
the test.
