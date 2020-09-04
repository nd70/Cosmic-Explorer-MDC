function params = readParamsFromFile(params)
%
%  readParamsFromFile --- read in search parameters from a file
%
%  readParamsFromFile(params) reads in search parameters from a
%  file, returning the parameters in a structure. The input structure
%  params may be partially filled in but must at least contain the 
%  field paramsFile specifying the full path to the params file.
%  If paramsFile is a null string it will be ignored and the params
%  structure will be returned unchanged.
%
%  Assumes parameters are given by name/value pair.
%
%  Routine written by Joseph D. Romano, John T. Whelan, Vuk Mandic.
%  Contact Joseph.Romano@astro.cf.ac.uk, john.whelan@ligo.org and/or
%  vmandic@ligo.caltech.edu
%
%  $Id: readParamsFromFile.m,v 1.38 2008-12-17 16:08:53 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% Check if the paramsFile is null or nonexistent
if (isempty(params.paramsFile))
  % Param file not specified, use defaults
  return;
elseif (~exist(params.paramsFile, 'file'))
  error(sprintf('Parameter file \''%s\'' not found.', params.paramsFile));
end;

%% read in name/value pairs
[names,values] = ...
  textread(params.paramsFile, '%s %s', -1, 'commentstyle', 'matlab', 'bufsize', 65536);

%% check that number of names and values are equal
if length(names) ~= length(values)
  error(sprintf('Invalid parameter file \''%s\'', number of variables does not match number of values.', params.paramsFile));
end

%% loop over parameter names, assigning values to structure
for ii=1:length(names)

  switch names{ii}

    case 'doDirectional'
      params.doDirectional = str2num(values{ii});

    case 'doSphericalHarmonics'
      params.doSphericalHarmonics = str2num(values{ii});

    case 'intermediate' %cet------------------------------
       params.intermediate = str2num(values{ii}); %-------

    case 'intFrameCachePath' %cet------------------------------
       params.intFrameCachePath = values{ii}; %-------

    case 'doNarrowbandRadiometer'
      params.doNarrowbandRadiometer = str2num(values{ii});

    case 'doAllSkyComparison'
      params.doAllSkyComparison = str2num(values{ii});

    case 'doFreqMask'
      params.doFreqMask = str2num(values{ii});

    case 'useDatafindServer'
      params.useDatafindServer = str2num(values{ii});

    case 'bkndstudy'
      params.bkndstudy = str2num(values{ii});

    case 'loudPixel'
      params.loudPixel = str2num(values{ii});

    case 'crunch_map'
      params.crunch_map = str2num(values{ii});

    case 'doBurstegard'
      params.doBurstegard = str2num(values{ii});

%% Begin backwards-compatibility option
    case 'doHighPass'
      params.doHighPass1 = str2num(values{ii});
      params.doHighPass2 = str2num(values{ii});
%% End backwards-compatibility option

    case 'doHighPass1'
      params.doHighPass1 = str2num(values{ii});

    case 'doHighPass2'
      params.doHighPass2 = str2num(values{ii});

    case 'doSidereal'
      params.doSidereal = str2num(values{ii});

    case 'doMonteCarlo'
      params.doMonteCarlo = str2num(values{ii});

    case 'doSimulatedPointSource'
      params.doSimulatedPointSource = str2num(values{ii});

    case 'doSimulatedSkyMap'
      params.doSimulatedSkyMap = str2num(values{ii});

    case 'doSimulatedDetectorNoise'
      params.doSimulatedDetectorNoise = str2num(values{ii});

    case 'doMCoffset'
      params.doMCoffset = str2num(values{ii});

    case 'doConstTimeShift'
      params.doConstTimeShift = str2num(values{ii});

    case 'doOverlap'
      params.doOverlap = str2num(values{ii});

    case 'doCombine'
      params.doCombine = str2num(values{ii});

    case 'doShift1'
      params.doShift1 = str2num(values{ii});

    case 'doShift2'
      params.doShift2 = str2num(values{ii});

    case 'doBadGPSTimes'
      params.doBadGPSTimes = str2num(values{ii});

    case 'heterodyned'
      params.heterodyned = str2num(values{ii});

    case 'writeResultsToScreen'
      params.writeResultsToScreen = str2num(values{ii});

    case 'writeStatsToFiles'
      params.writeStatsToFiles = str2num(values{ii});

    case 'writeOutputToMatFile'
      params.writeOutputToMatFile = str2num(values{ii});

%% "Naive" sigmas are those calculated from the PSD of the segment
%%  being analyzed, rather than with the sliding PSD estimator

    case 'writeNaiveSigmasToFiles'
      params.writeNaiveSigmasToFiles = str2num(values{ii});

    case 'writeSpectraToFiles'
      params.writeSpectraToFiles = str2num(values{ii});

    case 'writeSensIntsToFiles'
      params.writeSensIntsToFiles = str2num(values{ii});
      
    case 'writeCoherenceToFiles'
      params.writeCoherenceToFiles = str2num(values{ii});

    case 'writeCohFToFiles'
      params.writeCohFToFiles = str2num(values{ii});

    case 'writeOptimalFiltersToFiles'
      params.writeOptimalFiltersToFiles = str2num(values{ii});

    case 'writeOverlapReductionFunctionToFiles'
      params.writeOverlapReductionFunctionToFiles = str2num(values{ii});

    case 'writeCalPSD1sToFiles'
      params.writeCalPSD1sToFiles = str2num(values{ii});

    case 'writeCalPSD2sToFiles'
      params.writeCalPSD2sToFiles = str2num(values{ii});

    case 'ifo1'
      params.ifo1 = values{ii};

    case 'ifo2'
      params.ifo2 = values{ii};

    case 'doTimingTransientSubtraction1'
      params.doTimingTransientSubtraction1 = str2num(values{ii});

    case 'doTimingTransientSubtraction2'
      params.doTimingTransientSubtraction2 = str2num(values{ii});

    case 'TimingTransientFile1'
      params.TimingTransientFile1 = values{ii};

    case 'TimingTransientFile2'
      params.TimingTransientFile2 = values{ii};

    case 'site1'
      params.site1 = values{ii};

    case 'site2'
      params.site2 = values{ii};

    case 'azimuth1'
      params.azimuth1 = str2num(values{ii});

    case 'azimuth2'
      params.azimuth2 = str2num(values{ii});

    case 'segmentDuration'
      params.segmentDuration = str2num(values{ii});

    case 'numSegmentsPerInterval'
      params.numSegmentsPerInterval = str2num(values{ii});

    case 'ignoreMidSegment'
      params.ignoreMidSegment = str2num(values{ii});

    case 'deltaF'
      params.deltaF = str2num(values{ii});

    case 'flow'
      params.flow = str2num(values{ii});

    case 'fhigh'
      params.fhigh = str2num(values{ii});

    case 'alphaExp'
      params.alphaExp = str2num(values{ii});

    case 'fRef'
      params.fRef = str2num(values{ii});

    case 'maxSegmentsPerMatfile'
      params.maxSegmentsPerMatfile = str2num(values{ii});

    case 'useSignalSpectrumHfFromFile'
      params.useSignalSpectrumHfFromFile = str2num(values{ii});

    case 'SpHFreqIntFlag'
      params.SpHFreqIntFlag = str2num(values{ii});

    case 'IMScramblePhase'
      params.IMScramblePhase = str2num(values{ii});

    case 'HfFile'
      params.HfFile = values{ii};

    case 'HfFileInterpolateLogarithmic'
      params.HfFileInterpolateLogarithmic = str2num(values{ii});

    case 'gammaLM_coeffsPath'
      params.gammaLM_coeffsPath = values{ii};

    case 'SpHLmax'
      params.SpHLmax = str2num(values{ii});

    case 'useSkyPatternFile'
      params.useSkyPatternFile = str2num(values{ii});

    case 'SkyPatternFile'
      params.SkyPatternFile = values{ii};

    case 'SkyPatternRightAscensionNumPoints'
      params.SkyPatternRightAscensionNumPoints = str2num(values{ii});

    case 'SkyPatternDeclinationNumPoints'
      params.SkyPatternDeclinationNumPoints = str2num(values{ii});

    case 'maxCorrelationTimeShift'
      params.maxCorrelationTimeShift = str2num(values{ii});

    case 'UnphysicalTimeShift'
      params.UnphysicalTimeShift = str2num(values{ii});

    case 'minMCoff'
      params.minMCoff = str2num(values{ii});

    case 'maxMCoff'
      params.maxMCoff = str2num(values{ii});

    case 'ConstTimeShift'
      params.ConstTimeShift = str2num(values{ii});

    case 'ShiftTime1'
      params.ShiftTime1 = str2num(values{ii});

    case 'ShiftTime2'
      params.ShiftTime2 = str2num(values{ii});

    case 'maxDSigRatio'
      params.maxDSigRatio = str2num(values{ii});

    case 'minDSigRatio'
      params.minDSigRatio = str2num(values{ii});

    case 'minDataLoadLength'
      params.minDataLoadLength = str2num(values{ii});

%% Begin backwards-compatibility option
    case 'resampleRate'
      params.resampleRate1 = str2num(values{ii});
      params.resampleRate2 = str2num(values{ii});
%% End backwards-compatibility option

    case 'resampleRate1'
      params.resampleRate1 = str2num(values{ii});

    case 'resampleRate2'
      params.resampleRate2 = str2num(values{ii});

%% Begin KLUDGE because fshift not available in frames
    case 'fbase1'
      params.fbase1 = str2num(values{ii});

    case 'fbase2'
      params.fbase2 = str2num(values{ii});
%% End KLUDGE

%% Begin backwards-compatibility option
    case 'bufferSecs'
      params.bufferSecs1 = str2num(values{ii});
      params.bufferSecs2 = str2num(values{ii});
%% End backwards-compatibility option

    case 'bufferSecs1'
      params.bufferSecs1 = str2num(values{ii});

    case 'bufferSecs2'
      params.bufferSecs2 = str2num(values{ii});

%% Begin backwards-compatibility option
    case 'ASQchannel'
      params.ASQchannel1 = values{ii};
      params.ASQchannel2 = values{ii};
%% End backwards-compatibility option

    case 'ASQchannel1'
      params.ASQchannel1 = values{ii};

    case 'ASQchannel2'
      params.ASQchannel2 = values{ii};

    case 'injChannel1'
      params.injChannel1 = values{ii};

    case 'injChannel2'
      params.injChannel2 = values{ii};

    case 'injPrefix1'
      params.injPrefix1 = values{ii};

    case 'injPrefix2'
      params.injPrefix2 = values{ii};

%% Begin backwards-compatibility option
    case 'framePath'
      params.framePath1 = values{ii};
      params.framePath2 = values{ii};
%% End backwards-compatibility option

    case 'framePath1'
      params.framePath1 = values{ii};

    case 'framePath2'
      params.framePath2 = values{ii};

%% Begin backwards-compatibility option
    case 'frameType'
      params.frameType1 = values{ii};
      params.frameType2 = values{ii};
%% End backwards-compatibility option

    case 'frameType1'
      params.frameType1 = values{ii};

    case 'frameType2'
      params.frameType2 = values{ii};

%% Begin backwards-compatibility option
    case 'frameDuration'
      params.frameDuration1 = str2num(values{ii});
      params.frameDuration2 = str2num(values{ii});
%% End backwards-compatibility option

    case 'frameDuration1'
      params.frameDuration1 = str2num(values{ii});

    case 'frameDuration2'
      params.frameDuration2 = str2num(values{ii});

    case 'injFrameDuration1'
      params.injFrameDuration1 = str2num(values{ii});

    case 'injFrameDuration2'
      params.injFrameDuration2 = str2num(values{ii});

%% Begin backwards-compatibility option
    case 'hannDuration'
      params.hannDuration1 = str2num(values{ii});
      params.hannDuration2 = str2num(values{ii});
%% End backwards-compatibility option

    case 'hannDuration1'
      params.hannDuration1 = str2num(values{ii});

    case 'hannDuration2'
      params.hannDuration2 = str2num(values{ii});

%% Begin backwards-compatibility option
    case 'nResample'
      params.nResample1 = str2num(values{ii});
      params.nResample2 = str2num(values{ii});
%% End backwards-compatibility option

    case 'nResample1'
      params.nResample1 = str2num(values{ii});

    case 'nResample2'
      params.nResample2 = str2num(values{ii});

%% Begin backwards-compatibility option
    case 'betaParam'
      params.betaParam1 = str2num(values{ii});
      params.betaParam2 = str2num(values{ii});
%% End backwards-compatibility option

    case 'betaParam1'
      params.betaParam1 = str2num(values{ii});

    case 'betaParam2'
      params.betaParam2 = str2num(values{ii});

%% Begin backwards-compatibility option
    case 'highPassFreq'
      params.highPassFreq1 = str2num(values{ii});
      params.highPassFreq2 = str2num(values{ii});
%% End backwards-compatibility option

    case 'highPassFreq1'
      params.highPassFreq1 = str2num(values{ii});

    case 'highPassFreq2'
      params.highPassFreq2 = str2num(values{ii});

%% Begin backwards-compatibility option
    case 'highPassOrder'
      params.highPassOrder1 = str2num(values{ii});
      params.highPassOrder2 = str2num(values{ii});
%% End backwards-compatibility option

    case 'highPassOrder1'
      params.highPassOrder1 = str2num(values{ii});

    case 'highPassOrder2'
      params.highPassOrder2 = str2num(values{ii});

    % Flag to use cascade filter instead of direct filter
    case 'useCascadeFilter1'
      params.useCascadeFilter1 = str2num(values{ii});

    case 'useCascadeFilter2'
      params.useCascadeFilter2 = str2num(values{ii});

    case 'freqsToRemove'
      params.freqsToRemove = transpose(str2num(values{ii}));

    case 'nBinsToRemove'
      params.nBinsToRemove = transpose(str2num(values{ii}));

    case 'numTrials'
      params.numTrials = str2num(values{ii});

    case 'signalType'
      params.signalType = values{ii};

    case 'powerIndex'
      params.powerIndex = str2num(values{ii});

%% Begin backwards-compatibility option
    case 'simOmegaRef'
      params.simOmegaRef1 = str2num(values{ii});
      params.simOmegaRef2 = str2num(values{ii});
%% End backwards-compatibility option

    case 'simOmegaRef1'
      params.simOmegaRef1 = str2num(values{ii});

    case 'simOmegaRef2'
      params.simOmegaRef2 = str2num(values{ii});

    case 'simulationPath'
      params.simulationPath = values{ii};

    case 'simulatedPointSourcesFile'
      params.simulatedPointSourcesFile = values{ii};

    case 'simulatedPointSourcesPowerSpec'
      params.simulatedPointSourcesPowerSpec = values{ii};

    case 'simulatedPointSourcesInterpolateLogarithmic'
      params.simulatedPointSourcesInterpolateLogarithmic = str2num(values{ii});

    case 'simulatedPointSourcesBufferDepth'
      params.simulatedPointSourcesBufferDepth = str2num(values{ii});

    case 'simulatedPointSourcesHalfRefillLength'
      params.simulatedPointSourcesHalfRefillLength = str2num(values{ii});

    case 'simulatedPointSourcesNoRealData'
      params.simulatedPointSourcesNoRealData = str2num(values{ii});

    case 'simulatedPointSourcesMakeIncoherent'
      params.simulatedPointSourcesMakeIncoherent = str2num(values{ii});

    case 'simulatedSkyMapFile'
      params.simulatedSkyMapFile = values{ii};

    case 'simulatedSkyMapFileType'
      params.simulatedSkyMapFileType = str2num(values{ii});

    case 'simulatedSkyMapFileNumBins'
      params.simulatedSkyMapFileNumBins = str2num(values{ii});

    case 'simulatedSkyMapInjectAsSpH'
      params.simulatedSkyMapInjectAsSpH = str2num(values{ii});

    case 'simulatedSkyMapConvertLmax'
      params.simulatedSkyMapConvertLmax = str2num(values{ii});

    case 'simulatedSkyMapConvertDeg'
      params.simulatedSkyMapConvertDeg = str2num(values{ii});

    case 'simulatedSkyMapInjectTimeDomain'
      params.simulatedSkyMapInjectTimeDomain = str2num(values{ii});

    case 'simulatedSkyMapPowerSpec'
      params.simulatedSkyMapPowerSpec = values{ii};

    case 'simulatedSkyMapInterpolateLogarithmic'
      params.simulatedSkyMapInterpolateLogarithmic = str2num(values{ii});

    case 'simulatedSkyMapBufferDepth'
      params.simulatedSkyMapBufferDepth = str2num(values{ii});

    case 'simulatedSkyMapHalfRefillLength'
      params.simulatedSkyMapHalfRefillLength = str2num(values{ii});

    case 'simulatedSkyMapNoRealData'
      params.simulatedSkyMapNoRealData = str2num(values{ii});

    case 'simulatedDetectorNoisePowerSpec1'
      params.simulatedDetectorNoisePowerSpec1 = values{ii};

    case 'simulatedDetectorNoisePowerSpec2'
      params.simulatedDetectorNoisePowerSpec2 = values{ii};

    case 'simulatedDetectorNoiseInterpolateLogarithmic'
      params.simulatedDetectorNoiseInterpolateLogarithmic = str2num(values{ii});

    case 'simulatedDetectorNoiseBufferDepth'
      params.simulatedDetectorNoiseBufferDepth = str2num(values{ii});

    case 'simulatedDetectorNoiseHalfRefillLength'
      params.simulatedDetectorNoiseHalfRefillLength = str2num(values{ii});

    case 'simulatedDetectorNoiseNoRealData'
      params.simulatedDetectorNoiseNoRealData = str2num(values{ii});

    case 'badGPSTimesFile'
      params.badGPSTimesFile = values{ii};

    case 'alphaBetaFile1'
      params.alphaBetaFile1 = values{ii};

    case 'alphaBetaFile2'
      params.alphaBetaFile2 = values{ii};

    case 'calCavGainFile1'
      params.calCavGainFile1 = values{ii};

    case 'calCavGainFile2'
      params.calCavGainFile2 = values{ii};

    case 'calResponseFile1'
      params.calResponseFile1 = values{ii};

    case 'calResponseFile2'
      params.calResponseFile2 = values{ii};

    case 'gpsTimesPath1'
      params.gpsTimesPath1 = values{ii};

    case 'gpsTimesPath2'
      params.gpsTimesPath2 = values{ii};

    case 'frameCachePath1'
      params.frameCachePath1 = values{ii};

    case 'frameCachePath2'
      params.frameCachePath2 = values{ii};

    case 'doInjFromFile1'
      params.doInjFromFile1 = str2num(values{ii});

    case 'doInjFromFile2'
      params.doInjFromFile2 = str2num(values{ii});

    case 'injGPSTimesPath1'
      params.injGPSTimesPath1 = values{ii};

    case 'injGPSTimesPath2'
      params.injGPSTimesPath2 = values{ii};

    case 'injFrameCachePath1'
      params.injFrameCachePath1 = values{ii};

    case 'injFrameCachePath2'
      params.injFrameCachePath2 = values{ii};

    case 'injScale1'
      params.injScale1 = str2num(values{ii});

    case 'injScale2'
      params.injScale2 = str2num(values{ii});

    case 'outputDirPrefix'
      params.outputDirPrefix = values{ii};

    case 'outputFilePrefix'
      params.outputFilePrefix = values{ii};

    case 'suppressFrWarnings'
      params.suppressFrWarnings = str2num(values{ii});

    case 'jobsFileCommentStyle'
      params.jobsFileCommentStyle = values{ii};

    case 'verbose'
      params.verbose = str2num(values{ii});

    case 'saveGPS_diff'
      params.saveGPS_diff = str2num(values{ii});

    % stamp masking variables
    case 'doStampFreqMask'
     params.doStampFreqMask = str2num(values{ii});

    case 'StampFreqsToRemove'
     params.StampFreqsToRemove = str2num(values{ii});

    case 'StampnBinsToRemove'
     params.StampnBinsToRemove = str2num(values{ii});

    % eht on jan 5, 2011
    case 'DoStampInj'
      params.stamp_inj_doit = str2num(values{ii});

    case 'StampInjRA'
      params.stamp_inj_ra = str2num(values{ii});

    case 'StampInjDECL'
      params.stamp_inj_decl = str2num(values{ii});

    case 'StampInjStart'
      params.stamp_inj_start = str2num(values{ii});

    case 'StampInjType'
      params.stamp_inj_type = values{ii};

    case 'StampInjFile'
      params.stamp_inj_file = values{ii};

    case 'StampInjDur'
      params.stamp_inj_dur = str2num(values{ii});

    % used in preproc for STAMP studies-----------------
    case 'doDetectorNoiseSim'
      params.doDetectorNoiseSim = str2num(values{ii});

    case 'DetectorNoiseFile'
      params.DetectorNoiseFile = values{ii};
    %---------------------------------------------------

    % more params for stochmap addition
    case 'stochmap'
      params.stochmap = str2num(values{ii});

    case 'kludge'
      params.kludge = str2num(values{ii});

    case 'ra'
      params.ra = str2num(values{ii});

    case 'dec'
      params.dec = str2num(values{ii});

    case 'fmin'
      params.fmin = str2num(values{ii});

    case 'fmax'
      params.fmax = str2num(values{ii});

    case 'doPolar'
      params.doPolar = str2num(values{ii});

    case 'savePlots'
      params.savePlots = str2num(values{ii});

    case 'saveMat'
      params.saveMat = str2num(values{ii});

    case 'debug'
      params.debug = str2num(values{ii});

    case 'doMC'
      params.doMC = str2num(values{ii});

    case 'doRadon'
      params.doRadon = str2num(values{ii});

    case 'doBoxSearch'
      params.doBoxSearch = str2num(values{ii});

    case 'doClusterSearch'
      params.doClusterSearch = str2num(values{ii});

    case 'doLH'
      params.doLH = str2num(values{ii});

    case 'startGPS'
      params.startGPS = str2num(values{ii});

    case 'endGPS'
      params.endGPS = str2num(values{ii});

    case 'doRadiometer'
      params.doRadiometer = str2num(values{ii});

    case 'yMapScale'
      params.yMapScale = str2num(values{ii});

    case 'FMapScale'
      params.FMapScale = str2num(values{ii});

    case 'fluence'
      params.fluence = str2num(values{ii});

    case 'pp_seed'
      params.pp_seed = str2num(values{ii});

    case 'fft1dataWindow'
      params.fft1dataWindow = str2num(values{ii});
 
    case 'fft2dataWindow'
      params.fft2dataWindow = str2num(values{ii});

    case 'fixAntennaFactors'
      params.fixAntennaFactors = str2num(values{ii});

    case 'glitchCut'
      params.glitchCut = str2num(values{ii});

    case 'recordCuts'
      params.recordCuts = str2num(values{ii});

    case 'DQcut'
      params.DQcut = transpose(str2num(values{ii}));

    case 'phaseScramble'
      params.phaseScramble = str2num(values{ii});

    case 'purePlus'
      params.purePlus = str2num(values{ii});

    case 'psi'
      params.psi = str2num(values{ii});

    case 'iota'
      params.iota = str2num(values{ii});

    case 'bknd_study'
      params.bknd_study = str2num(values{ii});

    case 'bknd_study_dur'
      params.bknd_study_dur = str2num(values{ii});

    case 'pixelScramble'
      params.pixelScramble = str2num(values{ii});

    case 'Autopower'
      params.Autopower = str2num(values{ii});

    case 'powerinj'
      params.powerinj = str2num(values{ii});

    case 'injfile'
      params.injfile = str2num(values{ii});

    case 'inj_trials'
      params.inj_trials = str2num(values{ii});

    case 'alpha_n'
      params.alpha_n = str2num(values{ii});

    case 'alpha_max'
      params.alpha_max = str2num(values{ii});

    case 'alpha_min'
      params.alpha_min = str2num(values{ii});

    % Locust and Hough specific parameters
    case 'doLocust'
      params.doLocust = str2num(values{ii});

    case 'doHough'
      params.doHough = str2num(values{ii});

    case 'OutputDirectory'
      params.OutputDirectory = values{ii};

    case 'LocustCutThreshold'
      params.LocustCutThreshold = str2num(values{ii});

    case 'LocustRadiusX'
      params.LocustRadiusX = str2num(values{ii});

    case 'LocustRadiusY'
      params.LocustRadiusY = str2num(values{ii});

    case 'LocustRadSig'
      params.LocustRadSig = str2num(values{ii});

    case 'LocustDetectionThreshold'
      params.LocustDetectionThreshold = str2num(values{ii});

    case 'coarsegrained'
      params.coarsegrained = str2num(values{ii});

    case 'stamp_pem'
      params.stamp_pem = str2num(values{ii});

    case 'zeropad'
      params.zeropad = str2num(values{ii});

    case 'doDQT'
      params.doDQT = str2num(values{ii});

    case 'DQmatfile'
      params.DQmatfile = values{ii};

    case 'skypatch'
      params.skypatch = str2num(values{ii});

    case 'thetamax'
      params.thetamax = str2num(values{ii});

    case 'dtheta'
      params.dtheta = str2num(values{ii});

    case 'fastring'
      params.fastring = str2num(values{ii});

    case 'cluster.NN'
      params.cluster.NN = str2num(values{ii});

    case 'cluster.NR'
      params.cluster.NR = str2num(values{ii});

    case 'cluster.NCN'
      params.cluster.NCN = str2num(values{ii});

    case 'cluster.NCR'
      params.cluster.NCR = str2num(values{ii});

    case 'cluster.pixelThreshold'
      params.cluster.pixelThreshold = str2num(values{ii});

    case 'cluster.tmetric'
      params.cluster.tmetric = str2num(values{ii});

    case 'cluster.fmetric'
      params.cluster.fmetric = str2num(values{ii});
      
    case 'cluster.doCombineCluster'
      params.cluster.doCombineCluster = str2num(values{ii});

    case 'outputfilename'
      params.outputfilename = values{ii};

    case 'doCoincidentCut'
      params.doCoincidentCut = str2num(values{ii});

    case 'doMedianPSD'
      params.doMedianPSD = str2num(values{ii});


    otherwise
     warning(sprintf('Unknown parameter \''%s\'' in params file', names{ii}));

    end %% switch

end %% loop over parameter names

return
