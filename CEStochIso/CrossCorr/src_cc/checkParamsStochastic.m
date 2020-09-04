function checkParamsStochastic(params)
%
% Check parameters for consistency
% 
% input and output: params struct
%
% Routine copied from stochastic.m and modified my Stefan Ballmer 
% sballmer@caltech.edu 
%
% $Id: checkParamsStochastic.m,v 1.15 2009-05-20 13:31:02 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % Local variables for checking ranges
  maxFreq = 16384;
  maxGPSTime = 2147483647;

  if (params.startTime < 0 | params.startTime > maxGPSTime)
    error(sprintf('startTime %d is not a valid GPS time.', params.startTime));
  end;

  if (params.jobDuration < 0)
    error('jobDuration must be >= 0');
  end;

  if (params.startTime + params.jobDuration > maxGPSTime)
    error(sprintf('startTime + jobDuration %d is not a valid GPS time.', params.startTime + params.jobDuration));
  end;

  if (params.intermediate)
    checkPathExists('intFrameCachePath', params.intFrameCachePath);
    fprintf('Using intermediate data.\n');
  else
    if (~strcmp(params.gpsTimesPath1, '<auto>'))
      checkPathExists('gpsTimesPath1', params.gpsTimesPath1);
    end;

    if (~strcmp(params.gpsTimesPath2, '<auto>'))
      checkPathExists('gpsTimesPath2', params.gpsTimesPath2);
    end;

    checkPathExists('frameCachePath1', params.frameCachePath1);
    checkPathExists('frameCachePath2', params.frameCachePath2);
  end;

  if ~(strcmp(params.jobsFileCommentStyle,'matlab') ...
       | strcmp(params.jobsFileCommentStyle,'shell') ...
       | strcmp(params.jobsFileCommentStyle,'c') ...
       | strcmp(params.jobsFileCommentStyle,'c++'))
    error(sprintf('Unknown jobsFileCommentStyle \''%s\''.', params.jobsFileCommentStyle));
  end;

  if (params.doTimingTransientSubtraction1)
    checkFileExists('TimingTransientFile1', params.TimingTransientFile1);
  end;

  if (params.doTimingTransientSubtraction2)
    checkFileExists('TimingTransientFile2', params.TimingTransientFile2);
  end;

  if (params.segmentDuration <= 0)
    error(sprintf('segmentDuration (%.4e) is negative.', params.segmentDuration));
  end;

  if (params.numSegmentsPerInterval <= 0 | mod(params.numSegmentsPerInterval, 2) ~= 1)
    error(sprintf('numSegmentsPerInterval (%.4e) must be a positive odd number.', params.numSegmentsPerInterval));
  end;

  if (params.deltaF <= 0 | params.deltaF > maxFreq)
    error(sprintf('deltaF (%.4e) must be in the range (0, %.4e].', params.deltaF, maxFreq));
  end;

  % Possible for flow to be 0 (in testing for example)
  if (params.flow < 0 | params.flow > maxFreq)
    error(sprintf('flow (%.4e) must be in the range [0, %.4e].', params.flow, maxFreq));
  end;

  if (params.fhigh <= 0 | params.fhigh > maxFreq)
    error(sprintf('fhigh (%.4e) must be in the range (0, %.4e].', params.fhigh, maxFreq));
  end;

  if (params.fhigh < params.flow + params.deltaF)
    error(sprintf('fhigh (%.4e) is < flow + deltaF (%.4e).', ...
                  params.fhigh, params.flow + params.deltaF));
  end;

  if (params.fRef <= 0 | params.fRef > maxFreq)
    error(sprintf('fRef (%.4e) must be in the range (0, %.4e].', params.fRef, maxFreq));
  end;

  if params.useSignalSpectrumHfFromFile
    checkFileExists('HfFile,', params.HfFile);
  end;

  if (params.doSphericalHarmonics)
    checkPathExists('gammaLM_coeffsPath', params.gammaLM_coeffsPath);
  end;

  if (params.doDirectional)
    % no way to write the huge amount of data to screen
    if (params.writeResultsToScreen)
      warning(['doDirectional and writeResultsToScreen are both true, ', ...
               'this will produce a huge amount of output.' ]);
    end;

    if (params.useSkyPatternFile)
      checkFileExists('SkyPatternFile', params.SkyPatternFile);
    end;

    if (params.UnphysicalTimeShift ~= 0)
      warning('UnphysicalTimeShift is not 0, end result will be garbage.');
    end;
  end;

  if (params.useDatafindServer)
    server = getenv('LIGO_DATAFIND_SERVER');
    if (isempty(server))
      error('useDatafindServer is true but environment variable LIGO_DATAFIND_SERVER is not set.');
    end;

    [unixError, unixOut] = unix('/bin/sh -c "gw_data_find --ping 2>/dev/null"');
    if (unixError ~= 0)
      error(sprintf('\''gw_data_find --ping\'' fails for data grid server \''%s\'' specified in LIGO_DATAFIND_SERVER.', server));
    end;

   if (strcmp(params.frameType1, '<auto>'))
     error(sprintf('frameType1 may not be \''<auto>\'' when useDataFindServer is true.'));
   end;

   if (strcmp(params.frameType2, '<auto>'))
     error(sprintf('frameType2 may not be \''<auto>\'' when useDataFindServer is true.'));
   end;
  end;

  if (params.resampleRate1 <= 0 | params.resampleRate1 > maxFreq)
    error(sprintf('resampleRate1 (%.4e) must be in the range (0, %.4e].', params.resampleRate1, maxFreq));
  end;

  if (params.resampleRate2 <= 0 | params.resampleRate2 > maxFreq)
    error(sprintf('resampleRate2 (%.4e) must be in the range (0, %.4e].', params.resampleRate2, maxFreq));
  end;

  if (params.doInjFromFile1)
    % Not much else can be done here since we can't look in file to see
    % if the channel is present
    if (isempty(params.injChannel1))
      error(sprintf('injChannel1 is an empty string.'));
    end;

    if (~strcmp(params.injGPSTimesPath1, '<auto>'))
      checkFileExists('injGPSTimesPath1', params.injGPSTimesPath1);
    end;

    checkPathExists('injFrameCachePath1', params.injFrameCachePath1);
  end;

  if (params.doInjFromFile2)
    if (isempty(params.injChannel2))
      error(sprintf('injChannel2 is an empty string.'));
    end;

    if (~strcmp(params.injGPSTimesPath2, '<auto>'))
      checkFileExists('injGPSTimesPath2', params.injGPSTimesPath2);
    end;

    checkPathExists('injFrameCachePath2', params.injFrameCachePath2);
  end;

  if (params.hannDuration1 < 0 | params.hannDuration1 > params.segmentDuration)
    error(sprintf('hannDuration1 (%.4e) must be in the range (0, segmentDuration].', params.hannDuration1));
  end

  if (params.hannDuration2 < 0 | params.hannDuration2 > params.segmentDuration)
    error(sprintf('hannDuration2 (%.4e) must be in the range (0, segmentDuration].', params.hannDuration2));
  end

  if (params.highPassFreq1 <= 0 | params.highPassFreq1 > maxFreq)
    error(sprintf('highPassFreq1 (%.4e) must be in the range (0, %.4e].', params.highPassFreq1, maxFreq));
  end;

  if (params.highPassFreq2 <= 0 | params.highPassFreq2 > maxFreq)
    error(sprintf('highPassFreq2 (%.4e) must be in the range (0, %.4e].', params.highPassFreq2, maxFreq));
  end;

  if (params.doFreqMask)
    if (size(params.freqsToRemove, 2) ~= 1)
      error('freqsToRemove is not a column vector.');
    end;
    
    if (size(params.nBinsToRemove, 2) ~= 1)
      error('nBinsToRemove is not a column vector.');
    end;

    if (length(params.freqsToRemove) ~= length(params.nBinsToRemove))
      error(sprintf('Length of freqsToRemove (%d) is different to length of nBinsToRemove (%d).', ...
                    params.freqsToRemove, params.nBinsToRemove));
    end;
  end;
    
  if ~(strcmp(params.signalType,'const') | strcmp(params.signalType,'white'))
    error(sprintf('Unknown signalType \''%s\''.', params.signalType));
  end

  if params.doSimulatedPointSource
    checkPathExists('simulationPath', params.simulationPath);

    filename = [ params.simulationPath, '/', params.simulatedPointSourcesFile ];
    checkFileExists('simulatedPointSourcesFile', filename);

    filename = [ params.simulationPath, '/', params.simulatedPointSourcesPowerSpec ];
    checkFileExists('simulatedPointSourcesPowerSpec', filename);
  end;

  if params.doSimulatedSkyMap
    checkPathExists('simulationPath', params.simulationPath);
    
    filename = [ params.simulationPath, '/', params.simulatedSkyMapFile ];
    checkFileExists('simulatedSkyMapFile', filename);

    filename = [ params.simulationPath, '/', params.simulatedSkyMapPowerSpec ];
    checkFileExists('simulatedSkyMapPowerSpec', filename);

    clear filename;
  end;

  if params.doSimulatedDetectorNoise
    checkPathExists('simulationPath', params.simulationPath);

    filename = [ params.simulationPath, '/', params.simulatedDetectorNoisePowerSpec1 ];
    checkFileExists('simulatedDetectorNoisePowerSpec1', filename);

    filename = [ params.simulationPath, '/', params.simulatedDetectorNoisePowerSpec2 ];
    checkFileExists('simulatedDetectorNoisePowerSpec2', filename);

    clear filename;

    % is there an allowed range for minMCoff?
    if (params.maxMCoff < params.minMCoff)
      error(sprintf('maxMCoff (%.4e) is < minMCoff (%.4e)', params.maxMCoff, params.minMCoff));
    end;

    if (params.maxDSigRatio < params.minDSigRatio)
      error(sprintf('maxDSigRatio (%.4e) is < minDSigRatio (%.4e)', ...
                    params.maxDSigRatio, params.minDSigRatio));
    end;
  end;

  if (params.useCascadeFilter1)
    if (mod(params.highPassOrder1, 4) ~= 0)
      error(sprintf('highPassOrder1 (%d) must be divisible by 4 if cascade filter is used', params.highPassOrder1));
    end;
  else
    % A warning is emitted if a high-order filter is being used but not cascaded
    if (params.highPassOrder1 > 6)
      warning('checkParamsStochastic:filter', sprintf( 'highPassOrder1 (%d) > 6 but useCascadeFilter1 is false', params.highPassOrder1));
    end;
  end;

  if (params.useCascadeFilter2)
    if (mod(params.highPassOrder2, 4) ~= 0)
      error(sprintf('highPassOrder2 (%d) must be divisible by 4 if cascade filter is used', params.highPassOrder2));
    end;
  else
    % A warning is emitted if a high-order filter is being used but not cascaded
    if (params.highPassOrder2 > 6)
      warning('checkParamsStochastic:filter', sprintf( 'highPassOrder2 (%d) > 6 but useCascadeFilter2 is false', params.highPassOrder2));
    end;
  end;

  % It is fine for badGPSTimesFile to be a null string but if not we
  % need to check if the file exists
  if (~isempty(params.badGPSTimesFile))
    checkFileExists('badGPSTimesFile', params.badGPSTimesFile);
  end;

  % For these files 'none' is ok but anything else we must check
  if (~strcmp(params.alphaBetaFile1, 'none'))
    checkFileExists('alphaBetaFile1', params.alphaBetaFile1);
  end;

  if (~strcmp(params.alphaBetaFile2, 'none'))
    checkFileExists('alphaBetaFile2', params.alphaBetaFile2);
  end;

  if (~strcmp(params.calCavGainFile1, 'none'))
    checkFileExists('calCavGainFile1', params.calCavGainFile1);
  end;

  if (~strcmp(params.calCavGainFile2, 'none'))
    checkFileExists('calCavGainFile2', params.calCavGainFile2);
  end;

  if (~strcmp(params.calResponseFile1, 'none'))
    checkFileExists('calResponseFile1', params.calResponseFile1);
  end;

  if (~strcmp(params.calResponseFile2, 'none'))
    checkFileExists('calResponseFile2', params.calResponseFile2);
  end;

return;

