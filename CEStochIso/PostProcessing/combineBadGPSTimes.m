function combineBadGPSTimes(paramsFile, jobsFile, outputDir, badGPSTimesFile)
%
%  combineBadGPSTimes -- read the bad GPS times from the stochastic pipeline output
%                        and combine into a single file
%
% Inputs:
%   paramsFile - parameters file
%   jobsFile - list of science segments
%   outputDir - directory where post-processing results will go
%   badGPSTimesFile - Matlab file to which you wish to save badGPSTimes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % read in params structure from a file
  params = setStochasticParams(paramsFile, jobsFile, 1);
  checkParamsStochastic(params);

  % What if there are multiple trials?
  trial = 1;

  if (params.writeOutputToMatFile)
    badtimesfileprefix = [ params.outputFilePrefix '.job' ];
    fileSuffix = '.mat';
  else
    badtimesfileprefix = [ params.outputFilePrefix '_badGPSTimes.job' ];
    fileSuffix = [ '.trial' num2str(trial) '.dat' ];
  end;

  % Prefix for output from post-processing
  outputFileNamePrefix = [ outputDir '/' params.ifo1 params.ifo2 ];

  % read in job start time and duration from a file so we can count jobs and segments
  [startTimes, jobDurations] = readJobsFile(jobsFile, params.jobsFileCommentStyle);
  numJobs = length(startTimes);

  % Concatenate file names from all jobs into a cell array
  concatBadGPSTimes = {};
  for job = 1:numJobs
    concatBadGPSTimes{job} = [ badtimesfileprefix num2str(job) fileSuffix ];
  end;

  badGPSTimes = readBadGPSTimesFromFile(concatBadGPSTimes);

  save(badGPSTimesFile, 'badGPSTimes');

return;
