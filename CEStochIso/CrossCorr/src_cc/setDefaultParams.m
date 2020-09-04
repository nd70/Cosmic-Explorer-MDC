function params = setDefaultParams(paramsFile, startTime, jobDuration)
%
%  setDefaultParams --- set the default parameters, returning them in
%  a parameters structure.
%
%  Routine written by Joseph D. Romano, John T. Whelan, Vuk Mandic.
%  Contact Joseph.Romano@astro.cf.ac.uk, john.whelan@ligo.org and/or
%  vmandic@ligo.caltech.edu
%
%  $Id: $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % set the now time
  params.ddmmyyyyhhmmss = datestr(now);

  % add command line parameters to params struct
  params.paramsFile = paramsFile;

  % flags for optional operations
  params.doFreqMask = true;
  params.doHighPass1 = true;
  params.doHighPass2 = true;
  params.doMonteCarlo = false;
  params.doMCoffset = false;
  params.doConstTimeShift = false;
  params.doOverlap = true;
  params.heterodyned = false;
  params.suppressFrWarnings = false;
  params.doMedianPSD = false;

  % If true, the pipeline will attempt to source it's data from the LIGO data grid
  % server specified in the LIGO_DATAFIND_SERVER environment variable. In this
  % case the frame type must be specified and may not be '<auto>'.
  % If false, the pipeline will use the cache files found in frameCachePath.
  params.useDatafindServer = false;

  % ifo names
  params.ifo1 = 'H1';
  params.ifo2 = 'L1';

  % segment duration (sec)
  params.segmentDuration = 60;

  % parameters for sliding psd estimation:
  % numSegmentsPerInterval should be odd; ignoreMidSegment is a flag 
  % that allows you to ignore (if true) or include (if false) the 
  % analysis segment when estimating power spectra
  params.numSegmentsPerInterval = 3;
  params.ignoreMidSegment = true;

  % freq resolution and freq cutoffs for CC statistic sum (Hz)
  params.flow = 40;
  params.fhigh = 199.75;
  params.deltaF = 0.25;

  % params for Omega_gw (power-law exponent and reference freq in Hz)
  params.alphaExp = 0;
  params.fRef = 100;

  % resample rate (Hz)
  params.resampleRate1 = 1024;
  params.resampleRate2 = 1024;

  % buffer added to beginning and end of data segment to account for
  % filter transients (sec)
  params.bufferSecs1 = 2;
  params.bufferSecs2 = 2;

  % ASQ channel
  params.ASQchannel1 = 'LDAS-STRAIN';
  params.ASQchannel2 = 'LDAS-STRAIN';

  % frame type
  % When using gw_data_find to get frame data this must be set to the desired frame type.
  % It is ignored when using local cache files and may be set to '<auto>'.
  params.frameType1 = '<auto>';
  params.frameType2 = '<auto>';

  % frame duration in seconds. If duration is set to -1 the duration will be
  % deduced from the frame filename. Note that there is no proper validation
  % of frame duration because the frame reading is done by frgetvect which
  % does not return an error if the duration requested is longer than the
  % frame file duration, it just fills the rest with zeros. There doesn't appear
  % to be a way to determined the length of data in a given frame file using
  % the matlab routines provided in the frame library
  params.frameDuration1 = -1;
  params.frameDuration2 = -1;

  % duration of hann portion of tukey window 
  % (hannDuration = segmentDuration is a pure Hann window)
  params.hannDuration1 = 60;
  params.hannDuration2 = 60;

  % params for matlab resample routine
  params.nResample1 = 10;
  params.nResample2 = 10;
  params.betaParam1 = 5;
  params.betaParam2 = 5;

  % params for high-pass filtering (3db freq in Hz, and filter order) 
  params.highPassFreq1 = 32;
  params.highPassFreq2 = 32;
  params.highPassOrder1 = 6;
  params.highPassOrder2 = 6;

  % Default is to use direct filter rather than cascade
  params.useCascadeFilter1 = false;
  params.useCascadeFilter2 = false;

  % coherent freqs and number of freq bins to remove if doFreqMask=true;
  % NOTE: if an nBin=0, then no bins are removed even if doFreqMask=true
  % (coherent freqs are typically harmonics of the power line freq 60Hz
  % and the DAQ rate 16Hz)
  params.freqsToRemove = [71.9,80,96,108.3,109.9,112,120,128,144,150.75,152.25,160,176,180,192].';
  params.nBinsToRemove = [1,1,1,3,1,1,27,1,1,3,3,1,1,27,1].';

  % number of trials for monte carlo simulations (if doMonteCarlo = true)
  params.numTrials = 1;

  % type of SB signal to simulate 
  % (const for Omega_gw=const, white for Omega_gw propto f^3)
  params.signalType = 'const';

  % value of Omega_gw(f_Ref) for simulated SB signal
  params.simOmegaRef1 = 0;
  params.simOmegaRef2 = 0;

  % minimum and maximum time-shift for detector 1, if random offsets
  params.minMCoff = 0.1;
  params.maxMCoff = 0.6;

  % value of constant time-shift for detector 1, if used
  params.ConstTimeShift = 0.1;

  % calibration filenames
  params.alphaBetaFile1 = 'none';
  params.alphaBetaFile2 = 'none';
  params.calCavGainFile1 = 'none';
  params.calCavGainFile2 = 'none';
  params.calResponseFile1 = 'none';
  params.calResponseFile2 = 'none';

  % path to cache files
  params.gpsTimesPath1 = 'input/cachefiles/';
  params.gpsTimesPath2 = 'input/cachefiles/';
  params.frameCachePath1 = 'input/cachefiles/';
  params.frameCachePath2 = 'input/cachefiles/';

  % prefix for output filename
  params.outputFilePrefix = 'output/';

  %
  % Parameters moved from checkParamsStochastic
  %

  params.doSidereal = false;
  params.intermediate = false;
  params.intFrameCachePath = './intermediate/';
  params.jobsFileCommentStyle = 'matlab';

  % These parameters are consistent with doDirectional being false, they
  % will need to be set by the user in the params file if doDirectional
  % is set to true.
  % Note that by default writeResultsToScreen is true, which will produce a lot
  % of output if doDirectional is also true, so a warning has been added in
  % checkParamsStochastic for this situation
  params.doDirectional = false;
  params.useSkyPatternFile = false;
  params.SkyPatternFile = '';
  params.SkyPatternRightAscensionNumPoints = 0;
  params.SkyPatternDeclinationNumPoints = 0;
  params.maxCorrelationTimeShift = 0;
  params.UnphysicalTimeShift = 0;

  params.doSphericalHarmonics = false;

  % Must be set to a valid path if doSphericalHarmonics is true
  params.gammaLM_coeffsPath = '';
  params.doNarrowbandRadiometer = false;
  params.doAllSkyComparison = false;

  % These parameters are consistent with doSimulatedPointSource being false, they
  % will need to be set by the user in the params file if doSimulatedPointSource
  % is set to true
  params.doSimulatedPointSource = false;
  params.simulatedPointSourcesFile = '';
  params.simulatedPointSourcesPowerSpec = '';
  params.simulatedPointSourcesInterpolateLogarithmic = true;
  params.simulatedPointSourcesBufferDepth = 0;
  params.simulatedPointSourcesHalfRefillLength = 0;
  params.simulatedPointSourcesNoRealData = false;
  params.simulatedPointSourcesMakeIncoherent = 0;

  % This must be set to false to be consistent
  % with doDirectional being false by default  
  params.doAllSkyComparison = false;

  % These parameters are consistent with doSimulatedSkyMap being false.
  params.doSimulatedSkyMap = false;
  params.simulationPath = '';
  params.simulatedSkyMapFile = '';
  params.simulatedSkyMapFileType = 1;
  params.simulatedSkyMapFileNumBins = 1;
  params.simulatedSkyMapInjectAsSpH = true;
  params.simulatedSkyMapConvertLmax = 5;
  params.simulatedSkyMapConvertDeg = 10;
  params.simulatedSkyMapInjectTimeDomain = true;
  params.simulatedSkyMapPowerSpec = '';
  params.simulatedSkyMapInterpolateLogarithmic = true;
  params.simulatedSkyMapBufferDepth = 0;
  params.simulatedSkyMapHalfRefillLength = 0;
  params.simulatedSkyMapNoRealData = false;
  params.simulatedSkyMapMakeIncoherent = 0;

  % These parameters are consistent with doSimulatedDetectorNoise being false.
  params.doSimulatedDetectorNoise = false;
  params.simulatedDetectorNoisePowerSpec = '';  % for compatability
  params.simulatedDetectorNoisePowerSpec1 = '';
  params.simulatedDetectorNoisePowerSpec2 = '';
  params.simulatedDetectorNoiseInterpolateLogarithmic = true;
  params.simulatedDetectorNoiseMapBufferDepth = 0;
  params.simulatedDetectorNoiseMapHalfRefillLength = 0;
  % Change to true if doSimulatedDetectorNoise is true
  params.simulatedDetectorNoiseNoRealData = false;

  params.doMCoffset = false;
  params.doConstTimeShift = false;
  params.doOverlap = false;
  params.doCombine = false;
  params.doShift1 = false;
  params.doShift2 = false;
  params.doBadGPSTimes = false;

  params.doInjFromFile1 = false;
  % Must be set to a valid file if doInjFromFile1 is true
  params.injChannel1 = '';
  params.injFrameCachePath1 = '';
  params.injPrefix1 = params.ifo1;
  params.injFrameDuration1 = -1;
  params.injGPSTimesPath1 = '<auto>';
  params.injScale1 = 1;

  params.doInjFromFile2 = false;
  % Must be set to a valid file if doInjFromFile2 is true
  params.injChannel2 = '';
  params.injPrefix2 = params.ifo2;
  params.injFrameDuration2 = -1;
  params.injGPSTimesPath2 = '<auto>';
  params.injFrameCachePath2 = '';
  params.injScale2 = 1;

  params.writeResultsToScreen = true;
  params.writeStatsToFiles = true;
  params.writeOutputToMatFile = false;
  params.writeNaiveSigmasToFiles = false;
  params.writeSpectraToFiles = false;
  params.writeSensIntsToFiles = false;
  params.writeCoherenceToFiles = false;
  params.writeCohFToFiles = false;
  params.writeOptimalFiltersToFiles = false;
  params.writeOverlapReductionFunctionToFiles = false;
  params.writeCalPSD1sToFiles = false;
  params.writeCalPSD2sToFiles = false;
  params.doTimingTransientSubtraction1 = false;
  params.TimingTransientFile1 = '';
  params.doTimingTransientSubtraction2 = false;
  params.TimingTransientFile2 = '';
  params.azimuth1 = NaN;
  params.azimuth2 = NaN;
  params.maxSegmentsPerMatfile = 60;
  params.SpHFreqIntFlag = true;
  params.useSignalSpectrumHfFromFile = false;

  % Must be set to a valid filename if useSignalSpectrumHfFromFile is true
  params.HfFile = '';
  params.HfFileInterpolateLogarithmic = true;
  params.SpHLmax = 0;

  params.fbase1 = NaN;
  params.fbase2 = NaN;

  params.ShiftTime1 = 0;
  params.ShiftTime2 = 0;
  params.maxDSigRatio = 1000;
  params.minDSigRatio = 0;
  params.minDataLoadLength = 200;
  params.badGPSTimesFile = '';

  params.verbose = false;
  params.saveGPS_diff = false;

  % The third parameter to stochastic is interpreted as either the
  % job duration or the job number. Either way, it must be a number
  % or a numeric string.
  if (ischar(jobDuration) & length(str2num(jobDuration)) == 0)
    error('jobDuration must be numeric or a numeric string.');
  end;

  % Check which kind of parameters we have. If the second parameter
  % is a number or a string representing a number, treat it as the
  % start time, otherwise treat it as the name of the jobs file.
  if (isnumeric(startTime) | length(str2num(startTime)) ~= 0)
    % Set start time and duration based on parameters to stochastic.m
    params.startTime = strassign(startTime);
    params.jobDuration = strassign(jobDuration);
    params.jobsFile = '';
    params.jobNumber = params.startTime;
  else
    % Set start time and duration based on jobs file, interpreting
    % 'startTime' as jobsFile and 'jobDuration' as job number
    params.jobsFile = startTime;

    % convert string to numeric input argument (needed for compiled matlab) 
    params.jobNumber = strassign(jobDuration);

    % Read in job start time and duration from a file
    % This is one exceptional case that must be checked here
    % rather than in checkParamsStochastic
    checkFileExists('jobsFile', params.jobsFile);
    [startTimes, jobDurations] = readJobsFile(params.jobsFile, params.jobsFileCommentStyle);

    if (params.jobNumber < 0 | params.jobNumber > length(startTimes))
      error(sprintf('Job number %d outside range of available jobs %d', ...
                    params.jobNumber, length(startTimes)));
    else
      % if jobNumber is 0 it is a dummy job for testing
      if (params.jobNumber ~= 0)
        params.startTime = startTimes(params.jobNumber);
        params.jobDuration = jobDurations(params.jobNumber);
      else
        params.startTime = 0;
        params.jobDuration = 0;
      end;
    end;
  end;

  % if we need random numbers for each job
  randn('state', sum(100*clock)*params.jobNumber);

return;
