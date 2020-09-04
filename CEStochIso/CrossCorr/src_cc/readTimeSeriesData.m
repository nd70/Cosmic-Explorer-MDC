function [adcdata, dataOK] = readTimeSeriesData(channelName, dataStartTime, ...
                                                dataDuration, frameDuration, ...
                                                gpsTimesFile, frameCacheFile)
%
%  readTimeSeriesData --- reads in time-series data from frame files
%
%  readTimeSeriesData(channelName, dataStartTime, dataDuration, 
%                     frameDuration, gpsTimesFile, frameCacheFile)
%  uses a pre-fetched list of start times and file locations to find
%  and fetch the frame data.
%
%  If channelName starts with 'matlab:' the data are assumed to be in .mat
%  files rather than frames, and the string following the colon is assumed
%  to be the name of the variable containing the time series.
%
%  dataOK is a boolean (true or false) which indicates good or bad
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
  frameFiles = textread(frameCacheFile, '%s', -1, 'commentstyle', 'matlab');
  % parse the frame file names if necessary
  if (frameDuration <= 0 || strcmp(gpsTimesFile, '<auto>'))
    [gpsTimes, frameDurs] = decodeFrameNames(frameFiles);
  end;

  % if frame duration explicitly specified, use that instead
  % This overrides what is set by decodeFrameNames.
  if (frameDuration > 0)
    frameDurs = frameDuration*ones(length(frameFiles), 1);
  end;

  % if a GPS times file is specified instead of 'auto', read it. Note that
  % the gpsTimesFile has the form '<auto>gpsTimesX.N.txt' so we just examine
  % the first part of the filename.
  % This overrides what is set by decodeFrameNames.
  if (~strncmp(gpsTimesFile, '<auto>', length('<auto>')))
    gpsTimes = textread(gpsTimesFile, '%n', -1, 'commentstyle', 'matlab');
  end;

  magicLength = length('matlab:');
  trimmedName = channelName(4:end);
  if (strncmp(trimmedName, 'matlab:', magicLength))

    % get the data from Matlab files
    varName = trimmedName((magicLength+1):end);
    [adcdata, dataOK] = ...
        readTSDataFromMatfile(varName, dataStartTime, dataDuration, ...
                              gpsTimes, frameFiles, frameDurs);

  else

    % get the data from frame files
    chanObject = chanstruct(channelName);
    [vector, sampleRate, vectorError] = ...
        chanvector(chanObject, dataStartTime, dataDuration, gpsTimes, frameFiles, frameDurs);

    % check that the data is OK.
    % In order to make the pipeline more robust, we do not treat errors from chanvector as fatal to
    % the job at this time, although it is likely that some are.
    if (vectorError ~= 0)
      msg = translateChannelError(vectorError);
      warnMsg = [ 'Received error ''' msg ''' while trying to retrieve frame file list for' ...
                  ' channel ''' channelName  '''' ...
                  ' interval ' num2str(dataStartTime) '-' num2str(dataStartTime+dataDuration) '.' ];
      fprintf('READTIMESERIESDATA: %s\n', warnMsg);
      dataOK = false;
    end;

    if (dataOK )
      % fill time-series data structures
      adcdata.data   = transpose(vector);
      adcdata.tlow   = dataStartTime;
      adcdata.deltaT = 1/sampleRate;
    else
      % return all zeroes
      adcdata.data   = 0;
      adcdata.tlow   = 0;
      adcdata.deltaT = 0;
    end;

  end;

  % If the fbase or phase fields are not present, set them to
  % the dummy value NaN. This may be overridden in the params file.
  if (~isfield(adcdata,'fbase') )
    adcdata.fbase = NaN;
  end
  if (~isfield(adcdata,'phase') )
    adcdata.phase = NaN;
  end

return;
