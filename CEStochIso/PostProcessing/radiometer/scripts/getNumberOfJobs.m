[numJobs] = getNumberOfCombinedJobs(jobsToCombine, analysisFile)
%%%%%%%%%%%%%%%%%%%
% find maximum number of combined
% jobs from stochastic analysis.
% used in SumLocalRads.m
%%%%%%%%%%%%%%%%%%%

[paramsFileDir jobsFileDir] = textread(analysisFile, '%s %s',1);
[ifos epochs_vec skydir paramsFileName jobsFileName] = textread(analysisFile, '%s %s %s %s %s','headerlines',1);


