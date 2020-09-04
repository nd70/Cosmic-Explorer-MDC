function [adcdata, dataOK] = ...
   readTimeSeriesData_dev(channelName, dataStartTime, dataDuration, ...
		      frameType, frameDuration, gpsTimesFile, ...
			  frameCacheFile, NSegSkip, cache)
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
%  dev version, reads in different data, shifted by NSegSkip
%
%  $Id: readTimeSeriesData.m,v 1.9 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initially set dataOK to true
  dataOK = true;

% use a pre-fetched list of file locations to find and fetch the frame data
  old_frameFiles = textread(frameCacheFile, '%s\n', -1, ...
                       'commentstyle', 'matlab');
% parse frame filenames to get gpsTimes and frameDurs
  [old_gpsTimes,old_frameDurs]=decodeFrameNames(old_frameFiles);
  offset = dataStartTime-old_gpsTimes(1);

% EHT/
% assign new frameFiles based on new dataStartTime and dataDuration and using
% master cache
  index = [1:length(cache.frames2)];  % numerical index for frames2
% find the frame in cache.frames2 that matches the 1st frame in old_frameFiles
  old_first_frame = strcmp(cache.frames2, old_frameFiles(1))==1;
  % determine the index of that frame
  old_first_frame_index = sum(index'.*old_first_frame);
  NFrames = length(old_frameFiles);
  for ii=1:NFrames
    jj = old_first_frame_index+NSegSkip+(ii-1);
    % if new index exceeds array size, loop back around
    jj = mod(jj-1,length(cache.frames2))+1;
    frameFiles(ii) = cache.frames2(jj);
  end

  frameFiles = frameFiles';  % transpose array
  % parse frame filenames to get gpsTimes and frameDurs
  [gpsTimes,frameDurs]=decodeFrameNames(frameFiles);

  % now change dataStartTime
  dataStartTime = gpsTimes(1) + offset;
  % get the data
  chanObject  = chanstruct(channelName);
  [vector, sampleRate, vectorError] = chanvector(chanObject, ...
    dataStartTime, dataDuration, gpsTimes, frameFiles, frameDurs);

% EHT/
%  if necessary we can change the dataStartTime if there are missing data
%  issues though maybe this is not necessary
  if vectorError ~= 0
    keyboard
    dataStartTime = dataStartTime + 52;
  end
% /EHT 

  % check that the data is OK
  if vectorError == 0
    dataOK = true;
  else
    fprintf(['READTIMESERIESDATA: missing %s frame data starting' ...
             'at %d, ending at %d\n'], ...
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

return

