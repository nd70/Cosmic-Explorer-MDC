Getting Started with stochastic.m
by Eric Thrane (eric.thrane@ligo.org)

This file will walk you through the basic steps of running stochastic.m.  See 
also: https://wiki.ligo.org/Main/ComputingFAQ

-------------------------------------------------------------------------------
SETUP SCRIPTS AND MODULES
-------------------------------------------------------------------------------
* In order to run matlab successfully, make sure to include this command in
your .bash_login file as well as your condor executable wrapper(s):

     # ON CIT
     source FULL_PATH/cit-matlab-init

     # On Atlas
     source FULL_PATH/matlab_script_2012a.sh

* To run stochastic.m you will need to download part of matapps.  If you are
reading this file, you have already downloaded the CrossCorr module.  You will
also need some other matapps utilities.  I suggest that you use the ones that
we have set aside for STAMP:

   svn --username eric.thrane co https://ligo-vcs.phys.uwm.edu/svn/matapps/packages/stochastic/trunk/stamp2/

Downloading STAMP will automatically get you stochastic utilities such as 
frgetvect and coarsegrain.

* Next, you will need to set up your paths so that matlab knows where to find
the stochastic code suite.  To do this take a look at

    stochastic_paths.m

modify as necessary to point to your own directories, and then run every time
you open matlab.

* Before you open matlab, though, you will also need to run a setup script
from the UNIX shell:

  cit-matlab-init

-------------------------------------------------------------------------------
Test script
-------------------------------------------------------------------------------
With everything set up, we are now ready to run stochastic.m.

* stochastic.m can be demo'ed like this (from matlab):

  stochastic('example_params.txt', 'example_jobs.txt', 1);

* the first argument is the paramfile, which tells stochastic.m where to find
the raw data and how to process it.  Look at the example paramfile to see
what options you have

* the second argument is the jobfile.  the jobfile lists times associated with
uninterrupted data.  Look at the example jobfile to see more.

* The third argument is the job number.  In this tutorial, there are only three
jobs available.

* The paramfile points to cachefiles.  Cachefiles tell stochastic the paths to
the raw data.  You can see examples in cachefiles/.

* The output files of stochastic.m will all begin example_*.

-------------------------------------------------------------------------------
CREATING YOUR OWN JOBFILE AND CACHEFILES
-------------------------------------------------------------------------------
* If you are analyzing real data, try to get a jobfile from someone else who
has already analyzed the same data because they are a pain to make from
scratch.  If you have to make your own, you can use the jobfile scripts in 
../scripts/.

* If you are analyzing MC data, you can make your own jobfiles very easily by
writing a short script in perl or whatever your favorite language is.

* If you have a jobfile, but need new cachcefiles, this is fairly 
straightforward.  Just run the cachefiles.pl script in ../script/.  Before 
running any of these scripts, read the script to make sure that everything is
set correctly.

-------------------------------------------------------------------------------
COMPILING STOCHASTIC
-------------------------------------------------------------------------------
* Before you can condor submit stochastic jobs, you must compile stochastic.m.
Here's how.  First, from CrossCorr, create a directory called bin/ to build the
executable.  (This is for book-keeping purposes.)  Next, open matlab, and cd
to bin/.  Make sure all your matlab paths are corretly set (see above).  Then,
type

        mcc -R -nojvm -R -nodisplay -m stochastic

This will create an executable file in bin called stochastic (with no .m).  Try
running stochastic from the command line to make sure it runs.  You must run it
at least once from the command line anyway before condor submitting because the
first time a matlab executable runs, it must extract the ctf archive.

* Now you are ready to condor submit stochastic.  You will need

       -a wrapper (I call mine condor.sh) to source the matlab script before 
       calling stochastic).
       -a .dag file and/or a script to make your dag
       - a .sub file

These are all generic condor things, so we will not cover them in this readme!
However, there are condor resources on the sgwb wiki page:

         https://wiki.ligo.org/Main/StochasticGroup
