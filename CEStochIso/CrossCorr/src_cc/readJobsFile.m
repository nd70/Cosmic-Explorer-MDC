%
%  readJobsFile --- read the jobs file in a standardised way
%
%  Jobs file are generally in the format
%    <jobNumber> <GPS_start_time> <GPS_stop_time> <duration>
%
% with comments in Matlab format by default. This routine provides a
% standard method for reading in jobs files using the default comment
% style, which can be overridden by adding an optiona second parameter.
% The output is a vector containing the list of start times and another
% vector containing the list of durations. The job id and end time are
% not returned.
%
%  $Id: $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [startTimes, jobDurations]=readJobsFile(jobsFile, commentStyle)

  if (nargin < 2)
    commentStyle = 'matlab';
  end;

  [ignore1, startTimes, ignore2, jobDurations] = ...
      textread(jobsFile, '%n %n %n %n', -1, 'commentstyle', commentStyle);

return;