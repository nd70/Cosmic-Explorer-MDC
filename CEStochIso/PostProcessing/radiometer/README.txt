To run the radiometer post processing script in condor you need to specify two functions in your paramsfile:

1. exeDir (a directory specifying where you want your executable to be produced along with a run.sh script.)

2. dagFullName (full name of dag file and sub files to be produced, without their .dag or .sub extentions. preferrably specifying a full path but not necessarily.)

Otherwise, set doCondor = true and live the good life. 
