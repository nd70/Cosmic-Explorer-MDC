function [adcdata, dataOK] = ...
   readTimeSeriesData2(channelName, dataStartTime, dataDuration, ...
		      frameType, frameDuration, gpsTimesFile, ...
		      frameCacheFile, doDetectorNoiseSim)
%
%  readTimeSeriesData --- reads in time-series data from frame files
%
%  readTimeSeriesData(channelName, dataStartTime, dataDuration, 
%  frameType, frameDuration, gpsTimesFile, frameCacheFile)
%  uses a pre-fetched list of start times and file locations to find
%  and fetch the frame data.
%
%  If channelName starts with 'matlab:' the data are assumed to be in .mat
%  files rather than frames, and the string following the colon is assumed
%  to be the name of the variable containing the time series.
%
%  dataOK is a boolean (=true or false) which indicates good or bad
%  (i.e., missing) data
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org
%
%  $Id: readTimeSeriesData.m,v 1.9 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initially set dataOK to true
dataOK = true;

% use a pre-fetched list of file locations to find and fetch the
% frame data
frameFiles = textread(frameCacheFile, '%s\n', -1, ...
                       'commentstyle', 'matlab');

% parse the frame file names if necessary
if ((frameDuration <= 0) ...
    || (strncmp(gpsTimesFile,'<auto>',length('<auto>'))))
  [gpsTimes,frameDurs]=decodeFrameNames(frameFiles);
end;

% if frame duration explicitly specified, use that instead
if (frameDuration >= 0)
  frameDurs  = frameDuration*ones(length(frameFiles), 1);
end;

% if a GPS times file is specified instead of 'auto', read it
if (strncmp(gpsTimesFile,'<auto>',length('<auto>')) == false)
  gpsTimes   = textread(gpsTimesFile, '%n', -1, ...
                         'commentstyle', 'matlab');
end;

magicLength = length('matlab:');
trimmedName = channelName(4:end);
if (strncmp(trimmedName,'matlab:',magicLength))
  varName = trimmedName((magicLength+1):end);
  [adcdata, dataOK] = ...
      readTSDataFromMatfile(varName, dataStartTime, dataDuration, ...
			    gpsTimes, frameFiles, frameDurs);
else
  % get the data
  if (doDetectorNoiseSim)
    % Fill with simulated Gaussian noise.
     sampleRate = 16384;
     vector = randn(1,dataDuration*sampleRate);
     vectorError = 0;

  else
    % Get data from frames.
    chanObject  = chanstruct(channelName);
    [vector, sampleRate, vectorError] = chanvector(chanObject, ...
      dataStartTime, dataDuration, gpsTimes, frameFiles, frameDurs);
  end
  

  % check that the data is OK
  if vectorError == 0
    dataOK = true;
  else
    fprintf('READTIMESERIESDATA: missing %s frame data starting at %d, ending at %d\n', ...
            channelName(1), dataStartTime, dataStartTime+dataDuration);
    dataOK = false;
  end

  if dataOK 
    % fill time-series data structures
    adcdata.data   = transpose(vector);
    adcdata.tlow   = dataStartTime;
    adcdata.deltaT = 1/sampleRate;
  else
    % return all zeroes
    adcdata.data   = 0;
    adcdata.tlow   = 0;
    adcdata.deltaT = 0;
  end
end

return

