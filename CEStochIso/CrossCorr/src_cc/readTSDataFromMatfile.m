function [selecteddata, dataOK] = ...
      readTSDataFromMatfile(varName, dataStartTime, dataDuration,...
                            startTimes, matFiles, matDurs);
%
%  readTSDataFromMatfile --- read in time-series data from matfiles
%
%  readTSDataFromMatfile(dataStartTime, dataDuration, startTimes, 
%  matFiles, matDurs) reads in the time series data from the requested 
%  matfiles.
% 
%  Routine written by John T. Whelan.
%  Contact john.whelan@ligo.org
%
%  $Id: readTSDataFromMatfile.m,v 1.9 2006-10-25 16:46:26 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nFiles=length(startTimes);

if ( length(matFiles) ~= nFiles )
  error('number of GPS times does not match number of files');
end

if ( length(matDurs) ~= nFiles )
  error('number of GPS times does not match number of durations');
end

% Make sure the matlab data files are all consecutive and non-overlapping

for I=2:nFiles
  if ( startTimes(I) - startTimes(I-1) > matDurs(I-1) )
    error('Gap between matfiles');
  elseif ( startTimes(I) - startTimes(I-1) < matDurs(I-1) )
    error('Matfiles overlap');
  end  
end
endTimes = startTimes + matDurs;
dataEndTime = dataStartTime + dataDuration;

% Collect indices of all frame files relevant to data
includedFrameIndex = find((endTimes>dataStartTime)&(startTimes<dataEndTime));
if ( length(includedFrameIndex) < 1 )
  error('No matlab files to read\n');
end

% Read data and metadata from first matlab file
firstFrameIndex=includedFrameIndex(1);
load(matFiles{firstFrameIndex});
command = sprintf('selecteddata=%s;',varName);
eval(command);
sampleRate = length(selecteddata.data)/matDurs(firstFrameIndex);
%% fprintf('%d datapoints separated by %f = duration %f compared to expected %f\n',...
%% 	length(selecteddata.data), selecteddata.deltaT, ...
%% 	length(selecteddata.data)*selecteddata.deltaT, ...
%% 	matDurs(firstFrameIndex));
% TODO Should probably error check this (and lots of things about the metadata)

% Figure out which data point corresponds to start of series
offset = (dataStartTime - selecteddata.tlow) * sampleRate;
% trim unneeded start of series
selecteddata.data = selecteddata.data(1+offset:end);
selecteddata.tlow = dataStartTime;

% For backwards compatibility with v0 calibrated matfiles
if (isfield(selecteddata,'fhet') & ~ isfield(selecteddata,'fbase'))
  selecteddata.fbase = selecteddata.fhet;
  selecteddata = rmfield(selecteddata,'fhet');
end

if (isfield(selecteddata,'fbase') & ~ isnan(selecteddata.fbase))
  if ( selecteddata.fbase ~= floor(selecteddata.fbase) )
    error('non-integer heterodyne base frequency not implemented');
  end
end

% Concatenate data from second and remaining matlab files
% TODO: should check consistency of metadata
loopVec = includedFrameIndex(2:end);
% Make sure this is a row vector
loopVec = transpose(loopVec(:));
for I=loopVec
  load(matFiles{I});
  command = sprintf('adcdata=%s;',varName);
  eval(command);
  selecteddata.data = [selecteddata.data; adcdata.data];
end

% Trim off end of data
numdatapts = dataDuration * sampleRate;
selecteddata.data = selecteddata.data(1:numdatapts);

% Set this to true for now; this code is designed to throw an error
% rather than setting the flag to false.

dataOK = true;

return
