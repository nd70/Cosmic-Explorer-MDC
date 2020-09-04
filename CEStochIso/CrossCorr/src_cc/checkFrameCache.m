function checkFrameCache(params)
%
% checkFrameCache -- validates the frame cache and retrieves the cache from
%                    the specified LIGO data grid server if specified.
%
% The intent of this function is to validate the frame cache files and GPS times files
% before the pipeline commences.
%
%   - it checks if the frame cache file exists. If it doesn't but the user has set
%     useDatafindserver to true then it tries to retrieve it from LIGO_DATAFIND_SERVER
%   - if the GPS times pathis not '<auto>', it checks if the GPS times file exists.
%     If not it can retrieve it from LIGO_DATAFIND_SERVER as above.
%   - it verifies that the specified channel is present. Since this is an inefficient process
%     we only check the first file listed in the frame cache. It is unlikely that
%     the channel could be missing from other frame files but if it is that section of data
%     would be zero-filled later by readTimeSeriesData
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Only need to check this if we aren't using intermediate data
  if (~params.intermediate)
    checkFrameCacheFile(params.frameCacheFile1, ...
                        params.gpsTimesPath1, params.gpsTimesFile1, ...
                        params.ifo1(1), params.frameType1, params.channelName1, ...
                        params.startTime, params.jobDuration, params.useDatafindServer);
    
    checkFrameCacheFile(params.frameCacheFile2, ...
                        params.gpsTimesPath2, params.gpsTimesFile2, ...
                        params.ifo2(1), params.frameType2, params.channelName2, ...
                        params.startTime, params.jobDuration, params.useDatafindServer);
  end;

  % Also test injection frames if used
  if (params.doInjFromFile1)
    checkFrameCacheFile(params.injFrameCacheFile1, ...
                        params.injGPSTimesPath1, params.injFrameCacheFile1, ...
                        params.injPrefix1(1), '<auto>', params.injChannelName1, ...
                        params.startTime, params.jobDuration, false);
  end;

  if (params.doInjFromFile2)
    checkFrameCacheFile(params.injFrameCacheFile2, ...
                        params.injGPSTimesPath2, params.injFrameCacheFile2, ...
                        params.injPrefix2(1), '<auto>', params.injChannelName2, ...
                        params.startTime, params.jobDuration, false);
  end;

return;

function checkFrameCacheFile(frameCacheFile, gpsTimesPath, gpsTimesFile, site, frameType, channelName, startTime, jobDuration, useDatafindServer)

  cacheFileExists = exist(frameCacheFile, 'file');
  gpsTimesFileExists = exist(gpsTimesFile, 'file') | strcmp(gpsTimesPath, '<auto>');

  % If we are using a datafind server, we retrieve the frame list and GPS times if
  %  -- the frame cache file doesn't exist; or
  %  -- the GPS times file doesn't exist AND the GPS path is not '<auto>'
  %
  % If we aren't using a datafind server, it is an error if
  %  -- the cache file doesn't exist; or
  %  -- the GPS times file doesn't exist AND the GPS path is not '<auto>'
  if (useDatafindServer)
    
    % If either the cache file or GPS times file is absent, retrieve the frame list
    if (~cacheFileExists | ~gpsTimesFileExists)
      [ gpsTimes, frameFiles, listError, frameDurs ] = framelist(site, frameType, startTime, jobDuration);
      if (listError ~= 0)
        msg = translateChannelError(listError);
        errMsg = [ 'Received error ''' msg ''' while trying to retrieve frame file list for' ...
                   ' site ''' site '''' ...
                   ' frame type ''' frameType '''' ....
                   ' interval ' num2str(startTime) '-' num2str(startTime+jobDuration) '.' ];
        % When attempting to get the frame list we will consider any error to be fatal to the job
        error(errMsg);
      end;
    end;

    if (~cacheFileExists)
      createFrameCacheFile(frameCacheFile, frameFiles);
    end;

    if (~gpsTimesFileExists)
      createGPSTimesFile(gpsTimesFile, gpsTimes);
    end;

  else

    % We are not using a datafind server
    if (~cacheFileExists)
      error(sprintf('Frame cache file \''%s\'' not found.', frameCacheFile));
    end;
    if (~gpsTimesFileExists)
      error(sprintf('GPS times file \''%s\'' not found.', gpsTimesFile));
    end;

  end;

  % Check that the specified channel exists in the frame file. Since we have no
  % way to do this apart from attempting to read some data we will do it for
  % just the first frame file.
  frameFile = textread(frameCacheFile, '%s', 1, 'commentstyle', 'matlab');
  if (~channelExists(channelName, frameFile{1}))
    error(sprintf('Channel \''%s\'' does not exist in file \''%s\''.', ...
                  channelName, frameFile{1}));
  end;

return; % checkFrameCacheFile

function createFrameCacheFile(frameCacheFile, frameFiles)

  f = fopen(frameCacheFile, 'w');
  for k = 1:length(frameFiles)
    fprintf(f, '%s\n', frameFiles{k});
  end;
  fclose(f);

return;

function createGPSTimesFile(gpsTimesFile, gpsTimes)

  f = fopen(gpsTimesFile, 'w');
  for k = 1:length(gpsTimes)
    fprintf(f, '%d\n', gpsTimes(k));
  end;
  fclose(f);

return;
