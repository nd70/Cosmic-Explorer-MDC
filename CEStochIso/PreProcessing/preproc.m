function preproc_out = preproc(paramsFile, jobsFile, jobNumber)
%
%  preproc --- main routine for running the stochastic search
%
%  preproc(paramsFile, jobsFile, jobNumber) executes the stochastic
%  search specified by a set of parameters on the chosen job. The jobNumber
%  must be an integer from 0 up to the number of non-commented lines in
%  jobsFile. A jobNumber of 0 is a "dummy" job which can be used to
%  quickly run the executable when compiled under Matlab R14, and check that
%  the parameters are correct.
%
%  The parameter values are written to a file as are any or all of the
%  following:
%
%  - CC statistic values and theoretical sigmas
%  - naive theoretical sigmas (calculated from the analysis segment)
%  - CC spectra
%  - sensitivity integrands (the integrand of 1/theoretical variance)
%  - (complex) coherence between two channels
%  - optimal filter functions
%  - calibrated PSDs
%
%  Simulated stochastic background signals can also be injected into the
%  data if desired.
%
%  Routine written by Joseph D. Romano, John T. Whelan, Martin McHugh,
%  and Vuk Mandic.
%  Contact Joseph.Romano@astro.cf.ac.uk, john.whelan@ligo.org,
%  mmchugh@loyno.edu, and/or vmandic@ligo.caltech.edu
%
%  $Id: preproc.m,v 1.3 2007/06/01 00:41:00 vmandic Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tStart=tic;
ddmmyyyyhhmmss  = datestr(now);

% always define preproc_out to avoid compilation errors
preproc_out.done = false;

% convert string to numeric input argument (needed for compiled matlab) 
jobNumber = strassign(jobNumber);

% read in params structure from a file
params = readParamsFromFile(paramsFile);

% set random number generator.  use pp_seed<0 to randomize with clock.
try
  if params.pp_seed>0
    fprintf('setting seed to %i\n', params.pp_seed);
  else
    fprintf('setting seed using clock\n');
    params.pp_seed =  sum(100*clock)*jobNumber;
  end
catch
  fprintf('setting seed to 1\n');
  params.pp_seed = 200;
end
randn('state', params.pp_seed);

% #############################################################################
% extract parameters from the params structure setting to default values
% if necessary
try
  jobsFileCommentStyle = params.jobsFileCommentStyle;
catch
  jobsFileCommentStyle = 'matlab';
end;
try
  doSidereal=params.doSidereal;
catch
  doSidereal=false;
end
try
  doDirectional=params.doDirectional;
catch
  doDirectional=false;
end
try
  doNarrowbandRadiometer=params.doNarrowbandRadiometer;
catch
  doNarrowbandRadiometer=false;
end
try
  doAllSkyComparison=params.doAllSkyComparison;
catch
  doAllSkyComparison=false;
end
if ~doDirectional
  doAllSkyComparison=false;
end
try
  doFreqMask = params.doFreqMask;
catch
  error('doFreqMask parameter not set');
end
try
  doHighPass1 = params.doHighPass1;
catch
  error('doHighpass1 parameter not set');
end
try
  doHighPass2 = params.doHighPass2;
catch
 error('doHighpass2 parameter not set');
end
try
  doMonteCarlo = params.doMonteCarlo;
catch
  doMonteCarlo = false;
end
try
  doSimulatedPointSource = params.doSimulatedPointSource;
catch
  doSimulatedPointSource = false;
end
try
  doMCoffset = params.doMCoffset;
catch
  doMCoffset = false;
end
try
  doConstTimeShift = params.doConstTimeShift;
catch
  doConstTimeShift = false;
end
try
  doOverlap = params.doOverlap;
catch
  doOverlap = false;
end
try
  doCombine = params.doCombine;
catch
  doCombine = false;
end
try
  doShift1 = params.doShift1;
catch
  doShift1 = false;
end
try
  doShift2 = params.doShift2;
catch
  doShift2 = false;
end
try
  doBadGPSTimes = params.doBadGPSTimes;
catch
  doBadGPSTimes = false;
end
try
  heterodyned = params.heterodyned;
catch
  heterodyned = false;
end
try
  writeResultsToScreen = params.writeResultsToScreen;
catch
  writeResultsToScreen = false;
end
try
  writeStatsToFiles = params.writeStatsToFiles;
catch
  writeStatsToFiles = false;
end
try
  writeNaiveSigmasToFiles = params.writeNaiveSigmasToFiles;
catch
  writeNaiveSigmasToFiles = false;
end;
try
  writeSpectraToFiles = params.writeSpectraToFiles;
catch
  writeSpectraToFiles = false;
end
try
  writeSensIntsToFiles = params.writeSensIntsToFiles;
catch
  writeSensIntsToFiles = false;
end
try
    writeCoherenceToFiles = params.writeCoherenceToFiles;
catch
    writeCoherenceToFiles = false;
end
try
  writeOptimalFiltersToFiles = params.writeOptimalFiltersToFiles;
catch
  writeOptimalFiltersToFiles = false;
end
try
  writeOverlapReductionFunctionToFiles = params.writeOverlapReductionFunctionToFiles;
catch
  writeOverlapReductionFunctionToFiles = false;
end
%if not(doDirectional)
% currently sot supported (ORF is now unity...)
  writeOverlapReductionFunctionToFiles = false;
%end;
try
  writeCalPSD1sToFiles = params.writeCalPSD1sToFiles;
catch
  writeCalPSD1sToFiles = false;
end
try
  writeCalPSD2sToFiles = params.writeCalPSD2sToFiles;
catch
  writeCalPSD2sToFiles = false;
end
try
  ifo1 = params.ifo1;
catch
  error('ifo1 parameter not set');
end
try
  ifo2 = params.ifo2;
catch
  error('ifo2 parameter not set');
end
try
  doTimingTransientSubtraction1 = params.doTimingTransientSubtraction1;
  TimingTransientFile1          = params.TimingTransientFile1;
catch
  doTimingTransientSubtraction1 = false;
  TimingTransientFile1          = '';
end
try
  doTimingTransientSubtraction2 = params.doTimingTransientSubtraction2;
  TimingTransientFile2          = params.TimingTransientFile2;
catch
  doTimingTransientSubtraction2 = false;
  TimingTransientFile2          = '';
end
try
  azimuth1 = params.azimuth1;
catch
  azimuth1 = NaN;
end
try
  azimuth2 = params.azimuth2;
catch
  azimuth2 = NaN;
end
try
  segmentDuration = params.segmentDuration;
catch
  error('segmentDuration parameter not set');
end
try
  numSegmentsPerInterval = params.numSegmentsPerInterval;
catch
  error('numSegmentsPerInterval parameter not set');
end
try
  ignoreMidSegment = params.ignoreMidSegment;
catch
  error('ignoreMidSegment parameter not set');
end
try
  deltaF = params.deltaF;
catch
  error('deltaF parameter not set');
end
try
  flow = params.flow;
catch
  error('flow parameter not set');
end
try
  fhigh = params.fhigh;
catch
  error('fhigh parameter not set');
end
try
  alphaExp = params.alphaExp;
catch
  alphaExp = 0;
end
try
  fRef = params.fRef;
catch
  error('fRef parameter not set');
end
try
  maxSegmentsPerMatfile = params.maxSegmentsPerMatfile;
catch
  maxSegmentsPerMatfile = 60;
end
try
  useSignalSpectrumHfFromFile = params.useSignalSpectrumHfFromFile;
catch
  useSignalSpectrumHfFromFile = 0;
end
try
  HfFile = params.HfFile;
catch
  if useSignalSpectrumHfFromFile
    error('HfFile parameter not set');
  else
    HfFile = '';
  end
end
try
  HfFileInterpolateLogarithmic=params.HfFileInterpolateLogarithmic;
catch
  HfFileInterpolateLogarithmic=true;
end
if doDirectional
  % no way to write the uge amount of data to screen
  writeResultsToScreen=false;
  try useSkyPatternFile=params.useSkyPatternFile; catch useSkyPatternFile=false; end;
  try
    SkyPatternFile = params.SkyPatternFile;
  catch
    if useSkyPatternFile
      error('SkyPatternFile parameter not set');
    else
      SkyPatternFile = '';
    end
  end
  try
    SkyPatternRightAscensionNumPoints=params.SkyPatternRightAscensionNumPoints;
  catch
    error('SkyPatternRightAscensionNumPoints parameter not set');
  end;
  try
    SkyPatternDeclinationNumPoints=params.SkyPatternDeclinationNumPoints;
  catch
    error('SkyPatternDeclinationNumPoints parameter not set');
  end;
  try
    maxCorrelationTimeShift=params.maxCorrelationTimeShift;
  catch
    error('maxCorrelationTimeShift parameter not set');
  end;
  try
    UnphysicalTimeShift=params.UnphysicalTimeShift;
  catch
    UnphysicalTimeShift=0;
  end;
  if UnphysicalTimeShift~=0
    fprintf('Warning: UnphysicalTimeShift is not 0, end result will be garbage!\n');
  end
else
  useSkyPatternFile=0;
  SkyPatternFile='';
  SkyPatternRightAscensionNumPoints=0;
  SkyPatternDeclinationNumPoints=0;
  maxCorrelationTimeShift=0;
  UnphysicalTimeShift=0;
end;
try
  resampleRate1 = params.resampleRate1;
catch
  error('resampleRate1 parameter not set');
end
try
  resampleRate2 = params.resampleRate2;
catch
  error('resampleRate2 parameter not set');
end
try
  fbase1 = params.fbase1;
catch
  fbase1 = NaN;
end
try
  fbase2 = params.fbase2;
catch
  fbase2 = NaN;
end
try
  bufferSecs1 = params.bufferSecs1;
catch
  error('bufferSecs1 parameter not set');
end
try
  bufferSecs2 = params.bufferSecs2;
catch
  error('bufferSecs2 parameter not set');
end
try
  ASQchannel1 = params.ASQchannel1;
catch
  error('ASQchannel1 parameter not set');
end
try
  ASQchannel2 = params.ASQchannel2;
catch
  error('ASQchannel2 parameter not set');
end
try
  frameType1 = params.frameType1;
catch
  error('frameType1 parameter not set');
end
try
  frameType2 = params.frameType2;
catch
  error('frameType2 parameter not set');
end
try
  frameDuration1 = params.frameDuration1;
catch
  frameDuration1 = -1;
end
try
  frameDuration2 = params.frameDuration2;
catch
  frameDuration2 = -1;
end
try
  hannDuration1 = params.hannDuration1;
catch
  error('hannDuration1 parameter not set');
end
try
  hannDuration2 = params.hannDuration2;
catch
  error('hannDuration2 parameter not set');
end
try
  nResample1 = params.nResample1;
catch
  error('nResample1 parameter not set');
end
try
  nResample2 = params.nResample2;
catch
  error('nResample2 parameter not set');
end
try
  betaParam1 = params.betaParam1;
catch
  error('betaParam1 parameter not set');
end
try
  betaParam2 = params.betaParam2;
catch
  error('betaParam2 parameter not set');
end
try
  highPassFreq1 = params.highPassFreq1;
catch
  error('highPassFreq1 parameter not set');
end
try
  highPassFreq2 = params.highPassFreq2;
catch
  error('highPassFreq2 parameter not set');
end
try
  highPassOrder1 = params.highPassOrder1;
catch
  error('highPassOrder1 parameter not set');
end
try
  highPassOrder2 = params.highPassOrder2;
catch
  error('highPassOrder2 parameter not set');
end
try
  freqsToRemove = params.freqsToRemove;
catch
  error('freqsToRemove parameter not set');
end
try
  nBinsToRemove = params.nBinsToRemove;
catch
  error('nBinsToRemove parameter not set');
end
try
  numTrials = params.numTrials;
catch
  numTrials = 1;
end
try
  signalType = params.signalType;
catch
  signalType = 'const';
end
try
  simOmegaRef = params.simOmegaRef;
catch
  simOmegaRef = 0;
end
if doSimulatedPointSource
  try
    simulationPath=params.simulationPath;
  catch
    error('simulationPath parameter not set');
  end;
  try
    simulatedPointSourcesFile=params.simulatedPointSourcesFile;
  catch
    error('simulatedPointSourcesFile parameter not set');
  end;
  try
    simulatedPointSourcesPowerSpec=params.simulatedPointSourcesPowerSpec;
  catch
    error('simulatedPointSourcesPowerSpec parameter not set');
  end;
  try
    simulatedPointSourcesInterpolateLogarithmic=params.simulatedPointSourcesInterpolateLogarithmic;
  catch
    simulatedPointSourcesInterpolateLogarithmic=true;
  end;
  try
    simulatedPointSourcesBufferDepth=params.simulatedPointSourcesBufferDepth;
  catch
    error('simulatedPointSourcesBufferDepth parameter not set');
  end;
  try
    simulatedPointSourcesHalfRefillLength=params.simulatedPointSourcesHalfRefillLength;
  catch
    error('simulatedPointSourcesHalfRefillLength parameter not set');
  end;
  try
    simulatedPointSourcesNoRealData=params.simulatedPointSourcesNoRealData;
  catch
    simulatedPointSourcesNoRealData=false;
  end;
  try
    simulatedPointSourcesMakeIncoherent=params.simulatedPointSourcesMakeIncoherent;
  catch
    simulatedPointSourcesMakeIncoherent=0;
  end;
else
  simulationPath='';
  simulatedPointSourcesFile='';
  simulatedPointSourcesPowerSpec='';
  simulatedPointSourcesInterpolateLogarithmic=true;
  simulatedPointSourcesBufferDepth=0;
  simulatedPointSourcesHalfRefillLength=0;
  simulatedPointSourcesNoRealData=false;
  simulatedPointSourcesMakeIncoherent=0;
end;
try
  minMCoff = params.minMCoff;
catch
  minMCoff = 0.1;
end
try
  maxMCoff = params.maxMCoff;
catch
  maxMCoff = 0.6;
end
try
  ConstTimeShift = params.ConstTimeShift;
catch
  ConstTimeShift = 0;
end
try
  ShiftTime1 = params.ShiftTime1;
catch
  ShiftTime1 = 0;
end
try
  ShiftTime2 = params.ShiftTime2;
catch
  ShiftTime2 = 0;
end
try
  maxDSigRatio = params.maxDSigRatio;
catch
  maxDSigRatio = 1000;
end
try
  minDSigRatio = params.minDSigRatio;
catch
  minDSigRatio = 0;
end
try
  minDataLoadLength = params.minDataLoadLength;
catch
  minDataLoadLength = 200;
end
try
  badGPSTimesFile = params.badGPSTimesFile;
catch
  badGPSTimesFile = '';
end
try
  alphaBetaFile1 = params.alphaBetaFile1;
catch
  error('alphaBetaFile1 parameter not set');
end
try
  alphaBetaFile2 = params.alphaBetaFile2;
catch
  error('alphaBetaFile2 parameter not set');
end
try
  calCavGainFile1 = params.calCavGainFile1;
catch
  error('calCavGainFile1 parameter not set');
end
try
  calCavGainFile2 = params.calCavGainFile2;
catch
  error('calCavGainFile2 parameter not set');
end
try
  calResponseFile1 = params.calResponseFile1;
catch
  error('calResponseFile1 parameter not set');
end
try
  calResponseFile2 = params.calResponseFile2;
catch
  error('calResponseFile2 parameter not set');
end
try
  gpsTimesPath1 = params.gpsTimesPath1;
catch
  gpsTimesPath1 = '<auto>';
end
try
  gpsTimesPath2 = params.gpsTimesPath2;
catch
  gpsTimesPath2 = '<auto>';
end
try
  frameCachePath1 = params.frameCachePath1;
catch
  error('frameCachePath1 parameter not set');
end
try
  frameCachePath2 = params.frameCachePath2;
catch
  error('frameCachePath2 parameter not set');
end
try
  outputFilePrefix = params.outputFilePrefix;
catch
  error('outputFilePrefix parameter not set.\n');
end
try
  doDetectorNoiseSim = params.doDetectorNoiseSim;
catch
  doDetectorNoiseSim = 0;
end
try
  DetectorNoiseFile = params.DetectorNoiseFile;
catch
  DetectorNoiseFile = 'LIGOsrdPSD_40Hz.txt';
end
try
  doStochmap = params.stochmap;
  fprintf('doStochmap = %i\n', doStochmap);
catch
  doStochmap = 0;
end

try
  bknd_study = params.bknd_study;
catch
  bknd_study = 0;
end

if bknd_study~=0
  if ~exist('results/')
    error('You cannot run with bknd_study==true unless you have results/.');
  end
end

try
  bknd_study_dur = params.bknd_study_dur;
catch
  bknd_study_dur = 10;
end

try
  phaseScramble = params.phaseScramble;
catch
  phaseScramble = 0;
end

try
  pixelScramble = params.pixelScramble;
catch
  pixelScramble = 0;
end

% Define struct for STAMP injections.
try
  stamp_inj.doit = params.stamp_inj_doit;
  stamp_inj.ra = params.stamp_inj_ra;
  stamp_inj.decl = params.stamp_inj_decl;
  stamp_inj.start = params.stamp_inj_start;
  stamp_inj.type = params.stamp_inj_type;
  stamp_inj.file = params.stamp_inj_file;
catch
  stamp_inj.doit = 0;
end
if stamp_inj.doit
  switch stamp_inj.type
    case 'PSD'
      try
        stamp_inj.dur = params.stamp_inj_dur;
      catch
        error('StampInjDur must be set for injections of type PSD.\n');
      end
    case 'time_series'
    try
      stamp_inj.dur = params.stamp_inj_dur;
    catch
      stamp_inj.dur = -1;   % dur will be determined by injection file
    end
  end
end
% #############################################################################

% load injection data if requested.
if stamp_inj.doit
  fs = 16384;  % hard-coded LIGO sampling frequency
  dt = 1/fs;

  fprintf('loading STAMP injection "%s"...', stamp_inj.file);
  [h1 h2] = load_injection(params, stamp_inj, fs);
  fprintf('done.\n');

  % reset injection duration if it was changed by load_injection
  stamp_inj.dur = length(h1)/fs;

  % create array of times corresponding to injection time
  t_stamp0 = stamp_inj.start:1/fs:stamp_inj.start+stamp_inj.dur-1/fs;

  % pad the injection with zeros
  post_t = t_stamp0(end)+dt:dt:t_stamp0(end)+2*segmentDuration+2*bufferSecs1;
  postpad = zeros(size(post_t))';
  pre_t = t_stamp0(1)-2*segmentDuration-2*bufferSecs1 ...            
    : dt : t_stamp0(1)-dt;
  prepad = zeros(size(pre_t))';

  h1 = [prepad ; h1 ; postpad];
  h2 = [prepad ; h2 ; postpad];
  t_stamp = [pre_t t_stamp0 post_t];

  % The STAMP injection is shifted so as to compensate for the time shift.
  if params.doShift1
    StampShiftTime1 = -ShiftTime1;
  else
    StampShiftTime1 = 0;
  end
  if params.doShift2
    StampShiftTime2 = -ShiftTime2;
  else
    StampShiftTime2 = 0;
  end
  TShftMax = max(StampShiftTime1, StampShiftTime2);
  TShftMin = min(StampShiftTime1, StampShiftTime2);

end

% add few more things to the params structure
params.ddmmyyyyhhmmss = ddmmyyyyhhmmss;
params.paramsFile = paramsFile;
params.jobsFile = jobsFile;
params.matappsCVStag = 'matapps svn tag';

% determine bad GPS times, if given...
if doBadGPSTimes
  if isempty(badGPSTimesFile)
    badtimesstart = 9999999999;
    badtimesend = 0;
  else
%eht    [badtimesstart,badtimesend] = ctextread(badGPSTimesFile,'%f%f\n',-1,'commentstyle','matlab');
    [badtimesstart,badtimesend] = textread(badGPSTimesFile,'%f%f\n',-1,'commentstyle','matlab');
  end
end

% set total number of discrete frequencies
numFreqs = floor((fhigh-flow)/deltaF)+1;

% Take site letter (for detector geometry) from IFO name,
% unless overridden in parameter file
try
  site1 = params.site1;
catch
  site1 = getsitefromletter(ifo1(1));
end;
try
  site2 = params.site2;
catch
  site2 = getsitefromletter(ifo2(1));
end;
try
  suppressFrWarnings = params.suppressFrWarnings;
catch
  suppressFrWarnings = str2num('false');
end;
if suppressFrWarnings
  warning off frgetvect:info;
end;
try
  params.zeropad;
catch
  params.zeropad=true;
end;

% read in job start time and duration from a file
[ignore1, startTimes, ignore2, jobDurations] = ...
  textread(jobsFile, '%n %n %n %n', -1, ...
            'commentstyle', jobsFileCommentStyle);

% Check bounds of jobNumber. If jobNumber is zero,
% this is a dry run to verify parameters so it is
% not an error.
if (jobNumber == 0)
  warning('Dummy job number 0, all parameters verified - exiting');
  return;
end;
if (jobNumber < 0 | jobNumber > length(startTimes))
  error(sprintf('Job number %d outside range of available jobs %d', ...
    jobNumber, length(startTimes)));
end;

startTime   = startTimes(jobNumber);
jobDuration = jobDurations(jobNumber);

lastLoadedDataEnd1 = 0;
lastLoadedDataEnd2 = 0;

% get appropriate detector structure for each ifo
if (isnan(azimuth1))
  detector1 = getdetector(site1);
else
  detector1 = getdetector(site1,azimuth1);
end
if (isnan(azimuth2))
  detector2 = getdetector(site2);
else
  detector2 = getdetector(site2,azimuth2);
end

% construct filter coefficients for high-pass filtering
if doHighPass1
  [b1,a1] = butter(highPassOrder1, highPassFreq1/(resampleRate1/2), 'high');
end

if doHighPass2
  [b2,a2] = butter(highPassOrder2, highPassFreq2/(resampleRate2/2), 'high');
end

% set values for psd estimation (on resampled data, HP filtered data)
psdFFTLength1 = resampleRate1*(1/deltaF);
psdWindow1    = hann(psdFFTLength1);
psdOverlapLength1 = psdFFTLength1/2; 
detrendFlag1  = 'none';

psdFFTLength2 = resampleRate2*(1/deltaF);
psdWindow2    = hann(psdFFTLength2);
psdOverlapLength2 = psdFFTLength2/2; 
detrendFlag2  = 'none';

% set values for data windowing, zero-padding, and FFT
if doOverlap
  hannDuration1=segmentDuration;
  hannDuration2=segmentDuration;
end;
numPoints1    = segmentDuration*resampleRate1; 
dataWindow1   = tukeywin(numPoints1, hannDuration1/segmentDuration);
fftLength1    = 2*numPoints1;

numPoints2    = segmentDuration*resampleRate2;
dataWindow2   = tukeywin(numPoints2, hannDuration2/segmentDuration);
fftLength2    = 2*numPoints2;

% construct frequency mask for later use
data = constructFreqMask(flow, fhigh, deltaF, ...
                         freqsToRemove, nBinsToRemove, doFreqMask);
mask = constructFreqSeries(data, flow, deltaF);

% frequency vector and overlap reduction function
f = flow + deltaF*transpose([0:numFreqs-1]);
data  = overlapreductionfunction(f, detector1, detector2);
% If one of the data streams is heterodyned, set gamma.symmetry to 0
% so that negative frequencies are not included in the normalization.
if heterodyned
  gamma = constructFreqSeries(data, flow, deltaF, 0);
else
  gamma = constructFreqSeries(data, flow, deltaF, 1);
end

% construct name of gps times files and frame cache files for this job
gpsTimesFile1 = ...
  [gpsTimesPath1 'gpsTimes' ifo1(1) '.' num2str(jobNumber) '.txt'];
gpsTimesFile2 = ...
  [gpsTimesPath2 'gpsTimes' ifo2(1) '.' num2str(jobNumber) '.txt'];
frameCacheFile1 = ...
  [frameCachePath1 'frameFiles' ifo1(1) '.' num2str(jobNumber) '.txt'];
frameCacheFile2 = ...
  [frameCachePath2 'frameFiles' ifo2(1) '.' num2str(jobNumber) '.txt'];

% channel names
channelName1 = [ifo1 ':' ASQchannel1];
channelName2 = [ifo2 ':' ASQchannel2];

% read in calibration info

if ( ~strncmp(alphaBetaFile1,   'none', length(alphaBetaFile1))   & ...
     ~strncmp(calCavGainFile1,  'none', length(calCavGainFile1))  & ...
     ~strncmp(calResponseFile1, 'none', length(calResponseFile1)) )
[t1, f1, R01, C01, alpha1, gamma1] = ...
  readCalibrationFromFiles(alphaBetaFile1, calCavGainFile1, calResponseFile1);
end;

if ( ~strncmp(alphaBetaFile2,   'none', length(alphaBetaFile2))   & ...
     ~strncmp(calCavGainFile2,  'none', length(calCavGainFile2))  & ...
     ~strncmp(calResponseFile2, 'none', length(calResponseFile2)) )
[t2, f2, R02, C02, alpha2, gamma2] = ...
  readCalibrationFromFiles(alphaBetaFile2, calCavGainFile2, calResponseFile2);
end;

% check that numSegmentsPerInterval is odd    
if mod(numSegmentsPerInterval,2)==0
  error('numSegmentsPerInterval must be odd');
end

% determine number of intervals and segments to analyse
bufferSecsMax = max(bufferSecs1,bufferSecs2);


if ~doSidereal %do not worry about sidereal times
  M = floor( (jobDuration - 2*bufferSecsMax)/segmentDuration );

  if doOverlap
    numSegmentsTotal = 2*M-1;
    numIntervalsTotal = 2*(M - (numSegmentsPerInterval-1)) - 1;
    intervalTimeStride = segmentDuration/2;
  else 
    numSegmentsTotal = M;
    numIntervalsTotal = M - (numSegmentsPerInterval-1);
    intervalTimeStride = segmentDuration;
  end

  centeredStartTime = startTime + bufferSecsMax + ...
    floor( (jobDuration - 2*bufferSecsMax - M*segmentDuration)/ 2 );

else %calculate parameters compatible with the sidereal time
  srfac = 23.9344696 / 24; %sidereal time conversion factor
  srtime = GPStoGreenwichMeanSiderealTime(startTime) * 3600;
  md = mod(srtime,segmentDuration/srfac);
  centeredStartTime = round(startTime + segmentDuration - md*srfac); 

  %this is the start time of the first segment in this job that 
  %is compatible with the sidereal timing, but we still have to 
  %check that there is enough time for the preceding buffer

  if centeredStartTime - startTime < bufferSecsMax
    centeredStartTime = centeredStartTime + segmentDuration;
  end

  %now we can calculate the remaining bookkeeping variables
  M = floor( (jobDuration - bufferSecsMax - centeredStartTime + ...
              startTime) / segmentDuration );

  if params.doOverlap
    numIntervalsTotal = 2*(M - (numSegmentsPerInterval-1)) - 1;
    intervalTimeStride = segmentDuration/2;
  else 
    numIntervalsTotal = M - (numSegmentsPerInterval-1);
    intervalTimeStride = segmentDuration;
  end
end

 
% set number of trials to 1 if no monte carlo simulations
if doMonteCarlo==false; 
  numTrials = 1; 
end

isFirstPass=true;
% analyse the data
for I=1:numIntervalsTotal

  badSegmentData = false;
  badResponse = false;

  intervalStartTime = centeredStartTime + (I-1)*intervalTimeStride;

  % check if first pass through the loop
  if isFirstPass

    for J=1:numSegmentsPerInterval

      % read in time-series data from frames
      dataStartTime1 = intervalStartTime + (J-1)*segmentDuration - bufferSecs1;
      dataStartTime2 = intervalStartTime + (J-1)*segmentDuration - bufferSecs2;

      dataDuration1 = segmentDuration + 2*bufferSecs1;
      dataDuration2 = segmentDuration + 2*bufferSecs2;

      if dataStartTime1+dataDuration1 > lastLoadedDataEnd1
        lastLoadedDataEnd1 = min(dataStartTime1+minDataLoadLength,startTime+jobDuration);
        lastLoadedDataStart1 = dataStartTime1;
        tmpDuration = lastLoadedDataEnd1 - lastLoadedDataStart1;
        [longadcdata1, data1OK] = readTimeSeriesData2(channelName1,...
                              dataStartTime1, tmpDuration,...
                              frameType1, frameDuration1,...
			      gpsTimesFile1, frameCacheFile1,...
			      doDetectorNoiseSim);
      end


      if dataStartTime2+dataDuration2 > lastLoadedDataEnd2
        lastLoadedDataEnd2 = min(dataStartTime2+minDataLoadLength,startTime+jobDuration);
        lastLoadedDataStart2 = dataStartTime2;
        tmpDuration = lastLoadedDataEnd2 - lastLoadedDataStart2;
        [longadcdata2, data2OK] = readTimeSeriesData2(channelName2,...
                                dataStartTime2,tmpDuration,...
                                frameType2, frameDuration2,...
  		  	        gpsTimesFile2, frameCacheFile2,...
				doDetectorNoiseSim);
      end

      startindex = (dataStartTime1 - lastLoadedDataStart1) ...
                      /longadcdata1.deltaT + 1;
      endindex = (dataStartTime1 + dataDuration1 - lastLoadedDataStart1) ...
                      /longadcdata1.deltaT;
      adcdata1.data = longadcdata1.data(startindex:endindex);
      adcdata1.tlow = dataStartTime1;
      adcdata1.deltaT = longadcdata1.deltaT;

      startindex = (dataStartTime2 - lastLoadedDataStart2) ...
                      /longadcdata2.deltaT + 1;
      endindex = (dataStartTime2 + dataDuration2 - lastLoadedDataStart2) ...
                      /longadcdata2.deltaT;
      adcdata2.data = longadcdata2.data(startindex:endindex);
      adcdata2.tlow = dataStartTime2;
      adcdata2.deltaT = longadcdata2.deltaT;

      if doDetectorNoiseSim
         [adcdata1,no_need_1,no_need_2]=detectorNoise(adcdata1,DetectorNoiseFile);
         [adcdata2,no_need_1,no_need_2]=detectorNoise(adcdata2,DetectorNoiseFile);
         detectorNoise_data_1=adcdata1.data;
         detectorNoise_data_2=adcdata2.data;

%         adcdata1.data=zeros(size(adcdata1.data));
%         adcdata2.data=zeros(size(adcdata2.data));

      end

        % add STAMP injection if requested and if the time is appropriate
        if stamp_inj.doit
  	  if t_stamp0(end)+TShftMax>=adcdata1.tlow & ...
	    t_stamp0(1)+TShftMin < adcdata1.tlow+segmentDuration+2*bufferSecs1

             % crop the injection
            h1seg = h1(t_stamp+StampShiftTime1>=adcdata1.tlow & ...
              t_stamp+StampShiftTime1<adcdata1.tlow+segmentDuration+ ...
              2*bufferSecs1);
            h2seg = h2(t_stamp+StampShiftTime2>=adcdata2.tlow & ...
              t_stamp+StampShiftTime2<adcdata2.tlow+segmentDuration+ ...
              2*bufferSecs2);

            % add the injection
            adcdata1.data = adcdata1.data + h1seg;
            adcdata2.data = adcdata2.data + h2seg;

%      adcdata1.data = adcdata1.data + h1seg*1000;
%      adcdata2.data = adcdata2.data + h2seg*1000;

          end
        end

      % if either data stream is bad, set flag and exit loop
      if ( (data1OK==false) | (data2OK==false) )
        badSegmentData = true;
        break
      end

      if (~ isfield(adcdata1,'fbase') )
        adcdata1.fbase = NaN;
      end

      if (~ isfield(adcdata2,'fbase') )
        adcdata2.fbase = NaN;
      end

      if (~ isfield(adcdata1,'phase') )
	adcdata1.phase = NaN;
      end

      if (~ isfield(adcdata2,'phase') )
	adcdata2.phase = NaN;
      end

      % KLUDGE: can override base frequency in parameter file
      if (~ isnan(fbase1) )
	adcdata1.fbase = fbase1;
      end;
      if (~ isnan(fbase2) )
	adcdata2.fbase = fbase2;
      end;
      % End KLUDGE

      if ( isnan(adcdata1.fbase) & isnan(adcdata2.fbase) )
	if heterodyned
	  error('Trying to do heterodyned analysis on non-heterodyned data');
	end
      else
	if (~ heterodyned)
	  error('Trying to do non-heterodyned analysis on heterodyned data');
	end
      end

      % downsample the data 
      sampleRate1 = 1/adcdata1.deltaT;
      sampleRate2 = 1/adcdata2.deltaT;
      p1 = 1;  
      p2 = 1;  
      q1 = floor(sampleRate1/resampleRate1);
      q2 = floor(sampleRate2/resampleRate2);

      deltaT1 = 1/resampleRate1;
      deltaT2 = 1/resampleRate2;

      if sampleRate1 == resampleRate1
        data = adcdata1.data;
      else
        data = resample(adcdata1.data, p1, q1, nResample1, betaParam1);
      end
      n1(J) = constructTimeSeries(data, adcdata1.tlow, deltaT1, ...
                                  adcdata1.fbase, adcdata1.phase);

      if sampleRate2 == resampleRate2
        data = adcdata2.data;
      else
        data = resample(adcdata2.data, p2, q2, nResample2, betaParam2);
      end
      n2(J) = constructTimeSeries(data, adcdata2.tlow, deltaT2, ...
                                  adcdata2.fbase, adcdata2.phase);

      % free-up some memory
      clear adcdata1; 
      clear adcdata2;

      % calculate response functions from calibration data

      if ( strncmp(alphaBetaFile1,   'none', length(alphaBetaFile1))   | ...
           strncmp(calCavGainFile1,  'none', length(calCavGainFile1))  | ...
           strncmp(calResponseFile1, 'none', length(calResponseFile1)) )

        % the data is already calibrated
        response1(J) = constructFreqSeries(ones(numFreqs,1), flow, deltaF);
        transfer1(J) = response1(J);

      else
        calibsec1 = dataStartTime1 + bufferSecs1;
        [R1, responseOK1] = ...
          calculateResponse(t1, f1, R01, C01, alpha1, gamma1, calibsec1,...
			    ASQchannel1);

        % if response function is bad, set flag and exit loop
        if responseOK1==false
          badResponse = true;
          break
        end
  
        % evaluate response function at desired frequencies
        response1(J) = convertResponse(f1, R1, flow, deltaF, numFreqs, 0, 0);

        % convert to transfer functions (units: counts/strain) 
        transfer1(J) = convertResponse(f1, R1, flow, deltaF, numFreqs, 1, 0);
      end

      if ( strncmp(alphaBetaFile2,   'none', length(alphaBetaFile2))   | ...
           strncmp(calCavGainFile2,  'none', length(calCavGainFile2))  | ...
           strncmp(calResponseFile2, 'none', length(calResponseFile2)) )

        % the data is already calibrated
        response2(J) = constructFreqSeries(ones(numFreqs,1), flow, deltaF);
        transfer2(J) = response2(J);

      else
        calibsec2 = dataStartTime2 + bufferSecs2;
        [R2, responseOK2] = ...
          calculateResponse(t2, f2, R02, C02, alpha2, gamma2, calibsec2,...
			    ASQchannel2);

        % if response function is bad, set flag and exit loop
        if responseOK2==false
          badResponse = true;
          break
        end
   
        % evaluate response function at desired frequencies
        response2(J) = convertResponse(f2, R2, flow, deltaF, numFreqs, 0, 0);

        % convert to transfer functions (units: counts/strain) 
        transfer2(J) = convertResponse(f2, R2, flow, deltaF, numFreqs, 1, 0);
      end

    end % loop over segments J

    % if bad data or bad response function for any segment, continue with 
    % next interval
    if (badSegmentData | badResponse)
      continue
    else
      isFirstPass = false;
    end

  else

    % shift data and response functions accordingly
    for J=1:numSegmentsPerInterval-1

      if doOverlap
        % shift data by half a segment; need to worry about buffer
        N1 = length(n1(J).data);
        N2 = length(n2(J).data);
        bufferOffset1 = bufferSecs1/n1(J).deltaT;
        bufferOffset2 = bufferSecs2/n2(J).deltaT;

        data  = [n1(J).data(N1/2+1-bufferOffset1:N1-bufferOffset1) ; ...
                 n1(J+1).data(1+bufferOffset1:N1/2+bufferOffset1)];
        tlow  = n1(J).tlow+intervalTimeStride;
        n1(J) = constructTimeSeries(data, tlow, n1(J).deltaT, ...
                                    n1(J).fbase, n1(J).phase);

        data  = [n2(J).data(N2/2+1-bufferOffset2:N2-bufferOffset2) ; ...
                 n2(J+1).data(1+bufferOffset2:N2/2+bufferOffset2)];
        tlow  = n2(J).tlow+intervalTimeStride;
        n2(J) = constructTimeSeries(data, tlow, n2(J).deltaT, ...
                                    n2(J).fbase, n2(J).phase);

        % get response function corresponding to shifted start times

        if ( strncmp(alphaBetaFile1,   'none', length(alphaBetaFile1))   | ...
             strncmp(calCavGainFile1,  'none', length(calCavGainFile1))  | ...
             strncmp(calResponseFile1, 'none', length(calResponseFile1)) )

          % the data is already calibrated
          response1(J) = constructFreqSeries(ones(numFreqs,1), flow, deltaF);
          transfer1(J) = response1(J);

        else
          calibsec1 = n1(J).tlow + bufferSecs1;
          [R1, responseOK1] = ...
            calculateResponse(t1, f1, R01, C01, alpha1, gamma1, calibsec1,...
                            ASQchannel1);

          % if response function is bad, set flag and exit loop
          if responseOK1==false
            badResponse = true;
            break
          end
  
          % evaluate response function at desired frequencies
          response1(J) = convertResponse(f1, R1, flow, deltaF, numFreqs, 0, 0);

          % convert to transfer functions (units: counts/strain) 
          transfer1(J) = convertResponse(f1, R1, flow, deltaF, numFreqs, 1, 0);
        end

        if ( strncmp(alphaBetaFile2,   'none', length(alphaBetaFile2))   | ...
             strncmp(calCavGainFile2,  'none', length(calCavGainFile2))  | ...
             strncmp(calResponseFile2, 'none', length(calResponseFile2)) )

          % the data is already calibrated
          response2(J) = constructFreqSeries(ones(numFreqs,1), flow, deltaF);
          transfer2(J) = response2(J);

        else
          calibsec2 = n2(J).tlow + bufferSecs2;
          [R2, responseOK2] = ...
            calculateResponse(t2, f2, R02, C02, alpha2, gamma2, calibsec2,...
                            ASQchannel2);

          % if response function is bad, set flag and exit loop
          if responseOK2==false
            badResponse = true;
            break
          end
  
          % evaluate response function at desired frequencies
          response2(J) = convertResponse(f2, R2, flow, deltaF, numFreqs, 0, 0);

          % convert to transfer functions (units: counts/strain) 
          transfer2(J) = convertResponse(f2, R2, flow, deltaF, numFreqs, 1, 0);
        end

      else
        % simple shift by a full segment
        n1(J)=n1(J+1);
        n2(J)=n2(J+1);
     
        response1(J)=response1(J+1);
        response2(J)=response2(J+1);

        transfer1(J)=transfer1(J+1);
        transfer2(J)=transfer2(J+1);
      end

    end % loop over J

    % if bad response function for any segment, continue with next interval
    if badResponse
      continue
    end

    % read in time-series data for next segment
    dataStartTime1 = intervalStartTime ...
                     + (numSegmentsPerInterval-1)*segmentDuration ...
                     - bufferSecs1;
    dataStartTime2 = intervalStartTime ...
                     + (numSegmentsPerInterval-1)*segmentDuration ...
                     - bufferSecs2;

    dataDuration1  = segmentDuration + 2*bufferSecs1;
    dataDuration2  = segmentDuration + 2*bufferSecs2;

    if dataStartTime1+dataDuration1 > lastLoadedDataEnd1
      lastLoadedDataEnd1 = min(dataStartTime1+minDataLoadLength,startTime+jobDuration);
      lastLoadedDataStart1 = dataStartTime1;
      tmpDuration = lastLoadedDataEnd1 - lastLoadedDataStart1;
      [longadcdata1, data1OK] = readTimeSeriesData2(channelName1,...
                              dataStartTime1, tmpDuration,...
                              frameType1, frameDuration1,...
			      gpsTimesFile1, frameCacheFile1,...
			      doDetectorNoiseSim);
    end

    if dataStartTime2+dataDuration2 > lastLoadedDataEnd2
      lastLoadedDataEnd2 = min(dataStartTime2+minDataLoadLength,startTime+jobDuration);
      lastLoadedDataStart2 = dataStartTime2;
      tmpDuration = lastLoadedDataEnd2 - lastLoadedDataStart2;
      [longadcdata2, data2OK] = readTimeSeriesData2(channelName2,...
                                dataStartTime2,tmpDuration,...
                                frameType2, frameDuration2,...
	 		        gpsTimesFile2, frameCacheFile2,...
				doDetectorNoiseSim);
    end

    startindex = (dataStartTime1 - lastLoadedDataStart1) ...
                      /longadcdata1.deltaT + 1;
    endindex = (dataStartTime1 + dataDuration1 - lastLoadedDataStart1) ...
                      /longadcdata1.deltaT;
    adcdata1.data = longadcdata1.data(startindex:endindex);
    adcdata1.tlow = dataStartTime1;
    adcdata1.deltaT = longadcdata1.deltaT;

    startindex = (dataStartTime2 - lastLoadedDataStart2) ...
                      /longadcdata2.deltaT + 1;
    endindex = (dataStartTime2 + dataDuration2 - lastLoadedDataStart2) ...
                      /longadcdata2.deltaT;
    adcdata2.data = longadcdata2.data(startindex:endindex);
    adcdata2.tlow = dataStartTime2;
    adcdata2.deltaT = longadcdata2.deltaT;


      if doDetectorNoiseSim
         [adcdata1,no_need_1,no_need_2]=detectorNoise(adcdata1,DetectorNoiseFile);
         [adcdata2,no_need_1,no_need_2]=detectorNoise(adcdata2,DetectorNoiseFile);

         detectorNoise_data_1=adcdata1.data;
         detectorNoise_data_2=adcdata2.data;

         %find where the buffer ends (the overlap between old and new segment is where the first buffer ends through the end of the old segment)
         overlap_index_old_seg_1=bufferSecs1*sampleRate1+1;
         overlap_index_old_seg_2=bufferSecs2*sampleRate2+1;

         %find where the old data will be inserted in the new segment (beginning of segment up through buffer at end
         overlap_index_new_seg_1=(bufferSecs1+segmentDuration)*sampleRate1;
         overlap_index_new_seg_2=(bufferSecs2+segmentDuration)*sampleRate2; 

         adcdata1.data(1:overlap_index_new_seg_1)=detectorNoise_data_1(overlap_index_old_seg_1:end);
         adcdata2.data(1:overlap_index_new_seg_2)=detectorNoise_data_2(overlap_index_old_seg_2:end);         

      end

    % add STAMP injection if requested and if the time is appropriate
    if stamp_inj.doit
      if t_stamp0(end)+TShftMax>=adcdata1.tlow & ...
        t_stamp0(1)+TShftMin < adcdata1.tlow+segmentDuration+2*bufferSecs1

        h1seg = h1(t_stamp+StampShiftTime1>=adcdata1.tlow & ...
	  t_stamp+StampShiftTime1<adcdata1.tlow+segmentDuration+2*bufferSecs1);
        h2seg = h2(t_stamp+StampShiftTime2>=adcdata2.tlow & ...
	  t_stamp+StampShiftTime2<adcdata2.tlow+segmentDuration+2*bufferSecs2);

        % add the injection
        adcdata1.data = adcdata1.data + h1seg;
        adcdata2.data = adcdata2.data + h2seg;

 %     adcdata1.data = adcdata1.data + h1seg*1000;
 %     adcdata2.data = adcdata2.data + h2seg*1000;

      end
    end

    % if either data stream is bad, set flag and exit loop
    if ( (data1OK==false) | (data2OK==false) )
      badSegmentData = true;
      break
    end

    if (~ isfield(adcdata1,'fbase') )
      adcdata1.fbase = NaN;
    end
                                                                               
    if (~ isfield(adcdata2,'fbase') )
      adcdata2.fbase = NaN;
    end

    if (~ isfield(adcdata1,'phase') )
      adcdata1.phase = NaN;
    end

    if (~ isfield(adcdata2,'phase') )
      adcdata2.phase = NaN;
    end

    % KLUDGE: can override base frequency in parameter file
    if (~ isnan(fbase1) )
      adcdata1.fbase = fbase1;
    end
    if (~ isnan(fbase2) )
      adcdata2.fbase = fbase2;
    end
    % End KLUDGE
  
    if ( isnan(adcdata1.fbase) & isnan(adcdata2.fbase) )
      if heterodyned
        error('Trying to do heterodyned analysis on non-heterodyned data');
      end
    else
      if (~ heterodyned)
        error('Trying to do non-heterodyned analysis on heterodyned data');
      end
    end

   % downsample the data 
    sampleRate1 = 1/adcdata1.deltaT;
    sampleRate2 = 1/adcdata2.deltaT;
    p1 = 1;  
    p2 = 1;  
    q1 = floor(sampleRate1/resampleRate1);
    q2 = floor(sampleRate2/resampleRate2);

    deltaT1 = 1/resampleRate1;
    deltaT2 = 1/resampleRate2;

    if sampleRate1 == resampleRate1
      data = adcdata1.data;
    else
      data = resample(adcdata1.data, p1, q1, nResample1, betaParam1);
    end
    n1(numSegmentsPerInterval) = ...
      constructTimeSeries(data, adcdata1.tlow, deltaT1, ...
                          adcdata1.fbase, adcdata1.phase);

    if sampleRate2 == resampleRate2
      data = adcdata2.data;
    else
      data = resample(adcdata2.data, p2, q2, nResample2, betaParam2);
    end
    n2(numSegmentsPerInterval) = ...
      constructTimeSeries(data, adcdata2.tlow, deltaT2, ...
                          adcdata2.fbase, adcdata2.phase);

    % free-up some memory
    clear adcdata1; 
    clear adcdata2;

    % calculate response functions from calibration data

    if ( strncmp(alphaBetaFile1,   'none', length(alphaBetaFile1))   | ...
	 strncmp(calCavGainFile1,  'none', length(calCavGainFile1))  | ...
	 strncmp(calResponseFile1, 'none', length(calResponseFile1)) )

      % the data is already calibrated
      response1(numSegmentsPerInterval) = ...
        constructFreqSeries(ones(numFreqs,1), flow, deltaF);
      transfer1(numSegmentsPerInterval) = response1(numSegmentsPerInterval);

    else
      calibsec1 = dataStartTime1 + bufferSecs1;
      [R1, responseOK1] = ...
        calculateResponse(t1, f1, R01, C01, alpha1, gamma1, calibsec1,...
                            ASQchannel1);
    
    % if response function is bad, set flag and exit loop
    if responseOK1==false
      badResponse = true;
      break
    end
    %% evaluate response function at desired frequencies
    response1(numSegmentsPerInterval) = ...
      convertResponse(f1, R1, flow, deltaF, numFreqs, 0, 0);
    
      % convert to transfer functions (units: counts/strain) 
      transfer1(numSegmentsPerInterval) = ...
        convertResponse(f1, R1, flow, deltaF, numFreqs, 1, 0);
    end

    if ( strncmp(alphaBetaFile2,   'none', length(alphaBetaFile2))   | ...
         strncmp(calCavGainFile2,  'none', length(calCavGainFile2))  | ...
         strncmp(calResponseFile2, 'none', length(calResponseFile2)) )
	    
      % the data is already calibrated
      response2(numSegmentsPerInterval) = ...
	  constructFreqSeries(ones(numFreqs,1), flow, deltaF);
      transfer2(numSegmentsPerInterval) = response2(numSegmentsPerInterval);
      
    else
      calibsec2 = dataStartTime2 + bufferSecs2;
      [R2, responseOK2] = ...
        calculateResponse(t2, f2, R02, C02, alpha2, gamma2, calibsec2,...
                            ASQchannel2);
      
      % if response function is bad, set flag and exit loop
      if responseOK2==false
        badResponse = true;
        break
      end
  
      % evaluate response function at desired frequencies
      response2(numSegmentsPerInterval) = ...
        convertResponse(f2, R2, flow, deltaF, numFreqs, 0, 0);
    
      % convert to transfer functions (units: counts/strain) 
      transfer2(numSegmentsPerInterval) = ...
        convertResponse(f2, R2, flow, deltaF, numFreqs, 1, 0);
    end

  end % of if isFirstPass ... else ... end

  % if bad data or bad response function for any segment, continue with 
  % next interval
  if (badSegmentData | badResponse)
    continue
  end

  % initialize data array for average psds
  avg_data1 = zeros(numFreqs,1);
  avg_data2 = zeros(numFreqs,1);

  % loop over number of segments
  for J=1:numSegmentsPerInterval

    if doShift1
      shiftoffset = round(ShiftTime1 / n1(J).deltaT);
      qtempdata1 = circshift(n1(J).data,shiftoffset);
    else
      qtempdata1 = n1(J).data;
    end
    if doShift2
      shiftoffset = round(ShiftTime2 / n2(J).deltaT);
      qtempdata2 = circshift(n2(J).data,shiftoffset);
    else
      qtempdata2 = n2(J).data;
    end

    o1 = n1(J);
    o2 = n2(J);
    o1.data = qtempdata1;
    o2.data = qtempdata2;

    % high-pass filter the data (optional)
    if doHighPass1
      highpassed1 = constructTimeSeries(filtfilt(b1,a1,o1.data), ...
                                        o1.tlow, o1.deltaT, ...
		        		o1.fbase, o2.phase);
    else
      highpassed1 = o1;
    end
    if doHighPass2
      highpassed2 = constructTimeSeries(filtfilt(b2,a2,o2.data), ...
                                        o2.tlow, o2.deltaT, ...
					o2.fbase, o2.phase);
    else
      highpassed2 = o2;
    end

    % EHT on Jan 21, 2011: sometimes stamp injections create small imag parts
    highpassed1.data = real(highpassed1.data);
    highpassed2.data = real(highpassed2.data);

    % chop-off bad data at start and end of HP filtered, resampled data
    firstIndex1 = 1 + bufferSecs1*resampleRate1;
    firstIndex2 = 1 + bufferSecs2*resampleRate2;

    lastIndex1  = length(highpassed1.data)-bufferSecs1*resampleRate1;
    lastIndex2  = length(highpassed2.data)-bufferSecs2*resampleRate2;

    r1(J) = constructTimeSeries(highpassed1.data(firstIndex1:lastIndex1), ...
                                highpassed1.tlow + bufferSecs1, ...
                                highpassed1.deltaT, ...
				highpassed1.fbase, highpassed1.phase);
    r2(J) = constructTimeSeries(highpassed2.data(firstIndex2:lastIndex2), ...
                                highpassed2.tlow + bufferSecs2, ...
                                highpassed2.deltaT, ...
				highpassed2.fbase, highpassed2.phase);

    % estimate power spectra for optimal filter
    [temp1,freqs1] = psd(r1(J).data, psdFFTLength1, 1/r1(J).deltaT, ...
                         psdWindow1, psdOverlapLength1, detrendFlag1);
    [temp2,freqs2] = psd(r2(J).data, psdFFTLength2, 1/r2(J).deltaT, ...
                         psdWindow2, psdOverlapLength2, detrendFlag2);

    % normalize appropriately

    % If all the bins in the PSD are independent, we are dealing with
    % complex heterodyned data
    if length(temp1) == psdFFTLength1
      % Account for heterodyning of data.  This appears to be
      % the correct normalization because psd(), unlike pwelch(),
      % calculates the two-sided PSD of both real and complex data.
      freqs1shifted = fftshift(freqs1);
      spec1 = ...
          constructFreqSeries(2*r1(J).deltaT*fftshift(temp1), ...
		        	r1(J).fbase + freqs1shifted(1) ...
				- 1/r1(J).deltaT, ...
				freqs1(2)-freqs1(1), 0);
    else
      spec1 = ...
	  constructFreqSeries(2*r1(J).deltaT*temp1, freqs1(1), ...
				freqs1(2)-freqs1(1), 0);
    end

    % If all the bins in the PSD are independent, we are dealing with
    % complex heterodyned data
    if length(temp2) == psdFFTLength2
      % Account for heterodyning of data.  This appears to be
      % the correct normalization because psd(), unlike pwelch(),
      % calculates the two-sided PSD of both real and complex data.
      freqs2shifted = fftshift(freqs2);
      spec2 = ...
            constructFreqSeries(2*r2(J).deltaT*fftshift(temp2), ...
				r2(J).fbase + freqs2shifted(1) ...
				- 1/r2(J).deltaT, ...
				freqs2(2)-freqs2(1), 0);
    else
      spec2 = ...
	    constructFreqSeries(2*r2(J).deltaT*temp2, freqs2(1), ...
				freqs2(2)-freqs2(1), 0);
    end
  
    % coarse-grain noise power spectra to desired freqs
    psd1 = coarseGrain(spec1, flow, deltaF, numFreqs);
    psd2 = coarseGrain(spec2, flow, deltaF, numFreqs);

    % calibrate the power spectra
    calPSD1 = ...
      constructFreqSeries(psd1.data.*(abs(response1(J).data).^2), ...
                            psd1.flow, psd1.deltaF, psd1.symmetry);
    calPSD2 = ...
      constructFreqSeries(psd2.data.*(abs(response2(J).data).^2), ...
                            psd2.flow, psd2.deltaF, psd2.symmetry);
                        
    % calculate avg power spectra, ignoring middle segment if desired
    midSegment = (numSegmentsPerInterval+1)/2;
    if ( (ignoreMidSegment) & (J==midSegment) )
      % do nothing
      %fprintf('Ignoring middle segment\n');
    else
      avg_data1 = avg_data1 + calPSD1.data;
      avg_data2 = avg_data2 + calPSD2.data;
    end

    if (J==midSegment)
      % This calculates the "naive" theorerical variance, i.e.,
      % that calculated from the current segment without averaging 
      % over the whole interval.
      % This is useful for the stationarity veto which excludes
      % segments for which the naive sigma differs too much from
      % the one calculated with the sliding PSD average.
      [naiQ, naiVar, naiSensInt] ...
	    = calOptimalFilter(segmentDuration, gamma, ...
			       fRef, alphaExp, ... 
			       calPSD1, calPSD2, ...
			       dataWindow1, dataWindow2, mask);
      naiP1 = calPSD1.data;
      naiP2 = calPSD2.data;
    end
 
  end % loop over segments J

  % construct average power spectra
  if ignoreMidSegment
    avg_data1 = avg_data1/(numSegmentsPerInterval-1);
    avg_data2 = avg_data2/(numSegmentsPerInterval-1);
  else
    avg_data1 = avg_data1/numSegmentsPerInterval;
    avg_data2 = avg_data2/numSegmentsPerInterval;
  end

  calPSD1_avg = constructFreqSeries(avg_data1, flow, deltaF, 0);
  calPSD2_avg = constructFreqSeries(avg_data2, flow, deltaF, 0);

  % calculate optimal filter, theoretical variance, and sensitivity 
  % integrand using avg psds
  [Q, ccVar, sensInt] = calOptimalFilter(segmentDuration, gamma, ...
                                           fRef, alphaExp, ... 
                                           calPSD1_avg, calPSD2_avg, ...
                                           dataWindow1, dataWindow2, mask);

  % analyse the middle data segment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % window, zero-pad and fft the data
%  if params.zeropad
    rbartilde1 = windowAndFFT(r1(midSegment), dataWindow1, fftLength1);
    rbartilde2 = windowAndFFT(r2(midSegment), dataWindow2, fftLength2);
%  else
%  These lines are for diagnostic testing.
%    rbartilde1 = windowAndFFT_nopad(r1(midSegment), dataWindow1, fftLength1);
%    rbartilde2 = windowAndFFT_nopad(r2(midSegment), dataWindow2, fftLength2);
%  end

  % calculate the value and spectrum of the CC statistic
  [ccStat,CSD] = calCSD(rbartilde1, rbartilde2, Q, ...
                                     response1(midSegment), ...
                                     response2(midSegment));

  %recover the window parameters, add them to the structure
  [w1w2bar, w1w2squaredbar, w1w2ovlsquaredbar]=windowFactors(dataWindow1,...
                                                  dataWindow2);
  params.w1w2bar = w1w2bar;
  params.w1w2squaredbar = w1w2squaredbar;
  params.w1w2ovlsquaredbar = w1w2ovlsquaredbar;
  params.bias = 1/(2*(params.segmentDuration * params.deltaF * 2 - 1) * (9/11)) + 1;
  params.naivebias = 1/((params.segmentDuration * params.deltaF * 2 - 1) * (9/11)) + 1;
  params.SiderealTime = GPStoGreenwichMeanSiderealTime(r1(midSegment).tlow);
  

  %check if the interval is bad
  if doBadGPSTimes %determine if the current interval is bad
    c1 = intervalStartTime - bufferSecsMax < badtimesend;
    c2 = intervalStartTime + 3*segmentDuration + bufferSecsMax > badtimesstart;

    if sum(c1&c2) > 0 | sqrt(ccVar/naiVar)>maxDSigRatio | sqrt(ccVar/naiVar)<minDSigRatio
      badtimesflag = 1;
    else
      badtimesflag = 0;
    end
  else %do not check for bad times
    badtimesflag = 0;
  end

  %finally record the value in the relevant filename, if the bad time flag ok
  if badtimesflag == 0
    gpsstring = num2str(r1(midSegment).tlow);
%cet    filename = [outputDirPrefix outputFilePrefix gpsstring(1:4) '/' ...
    halfdur = segmentDuration/2;
%%cethrane Dec 6
%%    filename = [outputDirPrefix outputFilePrefix '-' ...
%%				gpsstring '-' num2str(halfdur) '.gwf'];
    filename = [outputFilePrefix '-' ...
				gpsstring '-' num2str(halfdur) '.gwf'];
    P1 = calPSD1_avg.data;
    P2 = calPSD2_avg.data;
    CC = 2 * CSD.data / w1w2bar / segmentDuration;
%    CC = CSD.data / w1w2bar;
%    save(filename,'CC','P1','P2','params');

    rec(1).channel = [ifo1 ':AdjacentPSD'];
    rec(1).data = P1;
    rec(1).type = 'd';
    rec(1).mode = 'a';

    rec(2).channel = [ifo2 ':AdjacentPSD'];
    rec(2).data = P2;
    rec(2).type = 'd';
    rec(2).mode = 'a';

    rec(3).channel = [ifo1 ':LocalPSD'];
    rec(3).data = naiP1;
    rec(3).type = 'd';
    rec(3).mode = 'a';

    rec(4).channel = [ifo2 ':LocalPSD'];
    rec(4).data = naiP2;
    rec(4).type = 'd';
    rec(4).mode = 'a';

%cet: change type = 'd' to 'dc' for complex data ... broken
%    rec(5).channel = [ifo1 ifo2 ':ReCSD'];
%    rec(5).data = real(CC);
%    rec(5).type = 'd';

    rec(5).channel = [ifo1 ifo2 ':CSD'];
    rec(5).data = CC;
    rec(5).type = 'dc';
    rec(5).mode = 'a';
%				    save('test.mat','CC'); %cet-Nov18 test
%cet-------try to pass parameters as such---------------
%    etparams.name='ParametersChannelName';
    etparams.name = [ifo1 ifo2 ':Params'];
%-------------------------------------------------------    
    if ~(doStochmap)
      doublestr = [];
      tt = fieldnames(params);
      kk=11; %cet
      for jj = 1:length(tt)
        tmp = getfield(params,tt{jj});
        if ischar(tmp)
          tmpstr = [tt{jj} ' ' tmp];
        else
          if size(tmp,1)>1
            %if it is an array, needs to be a row 
            tmp = transpose(tmp);
          end
          tmpstr = [tt{jj} ' ' num2str(tmp)];
          tmp = num2str(tmp); %cet: 10-10
        end
%cet-------pass parameters using new tool---------------
%%cethrane Dec 6
%%      etparams.parameters{jj*2-1}=tt{jj};
%%      etparams.parameters{jj*2}=tmp;
        ttvals{jj}=tmp;
% note - values that are equal to false are changed to 0
%-------------------------------------------------------
        doublestr = [doublestr double(tmpstr) 10];
      end
      %%cethrane Dec 6
      etparams.parameters=horzcat(tt, ttvals');
    end

%    rec(7).channel = [ifo1 ifo2 ':Params'];
%    rec(7).data = doublestr;
%    rec(7).type = 'd';
%    rec(7).mode = 'a';

%cet: all subsequent rec's shifted down by 1

    rec(6).channel = [ifo1 ifo2 ':flow'];
    rec(6).data = params.flow;
    rec(6).type = 'd';
    rec(6).mode = 'a';

    rec(7).channel = [ifo1 ifo2 ':fhigh'];
    rec(7).data = params.fhigh;
    rec(7).type = 'd';
    rec(7).mode = 'a';

    rec(8).channel = [ifo1 ifo2 ':deltaF'];
    rec(8).data = params.deltaF;
    rec(8).type = 'd';
    rec(8).mode = 'a';

    rec(9).channel = [ifo1 ifo2 ':GPStime'];
    rec(9).data = r1(midSegment).tlow;
    rec(9).type = 'd';
    rec(9).mode = 'a';

    rec(10).channel = [ifo1 ifo2 ':segmentDuration'];
    rec(10).data = params.segmentDuration;
    rec(10).type = 'd';
    rec(10).mode = 'a';

    rec(11).channel = [ifo1 ifo2 ':w1w2bar'];
    rec(11).data = params.w1w2bar;
    rec(11).type = 'd';
    rec(11).mode = 'a';

    rec(12).channel = [ifo1 ifo2 ':w1w2squaredbar'];
    rec(12).data = params.w1w2squaredbar;
    rec(12).type = 'd';
    rec(12).mode = 'a';
				    
    rec(13).channel = [ifo1 ifo2 ':w1w2ovlsquaredbar'];
    rec(13).data = params.w1w2ovlsquaredbar;
    rec(13).type = 'd';
    rec(13).mode = 'a';

    rec(14).channel = [ifo1 ifo2 ':bias'];
    rec(14).data = params.bias;
    rec(14).type = 'd';
    rec(14).mode = 'a';
 
    rec(15).channel = [ifo1 ifo2 ':naivebias'];
    rec(15).data = params.naivebias;
    rec(15).type = 'd';
    rec(15).mode = 'a';
 
    rec(16).channel = [ifo1 ifo2 ':SiderealTime'];
    rec(16).data = params.SiderealTime;
    rec(16).type = 'd';
    rec(16).mode = 'a';

    rec(17).channel = [ifo1 ifo2 ':numSegmentsPerInterval'];
    rec(17).data = params.numSegmentsPerInterval;
    rec(17).type = 'd';
    rec(17).mode = 'a';

    rec(18).channel = [ifo1 ifo2 ':resampleRate1'];
    rec(18).data = resampleRate1;
    rec(18).type = 'd';
    rec(18).mode = 'a';

    rec(19).channel = [ifo1 ifo2 ':resampleRate2'];
    rec(19).data = resampleRate2;
    rec(19).type = 'd';
    rec(19).mode = 'a';

    rec(20).channel = [ifo1 ifo2 ':nResample1'];
    rec(20).data = nResample1;
    rec(20).type = 'd';
    rec(20).mode = 'a';

    rec(21).channel = [ifo1 ifo2 ':nResample2'];
    rec(21).data = nResample2;
    rec(21).type = 'd';
    rec(21).mode = 'a';

    rec(22).channel = [ifo1 ifo2 ':betaParam1'];
    rec(22).data = betaParam1;
    rec(22).type = 'd';
    rec(22).mode = 'a';

    rec(23).channel = [ifo1 ifo2 ':betaParam2'];
    rec(23).data = betaParam2;
    rec(23).type = 'd';
    rec(23).mode = 'a';

    if ~(doStochmap)
      mkframe(filename,rec,'n',segmentDuration,r1(midSegment).tlow,etparams);
    else
      params.pass.P1(:,I) = rec(1).data;
      params.pass.P2(:,I) = rec(2).data;
      params.pass.naiP1(:,I) = rec(3).data;
      params.pass.naiP2(:,I) = rec(4).data;
      params.pass.CC(:,I) = rec(5).data;
      params.pass.ppflow(:,I) = rec(6).data;
      params.pass.ppfhigh(:,I) = rec(7).data;
      params.pass.ppdeltaF(:,I) = rec(8).data;
      params.pass.GPStimes(:,I) = rec(9).data;
      params.pass.segmentDuration(:,I) = rec(10).data;
      params.pass.ppw1w2bar(:,I) = rec(11).data;
      params.pass.ppw1w2squaredbar(:,I) = rec(12).data;
      params.pass.ppw1w2ovlsquaredbar(:,I) = rec(13).data;
      params.pass.bias(:,I) = rec(14).data;
      params.pass.naivebias(:,I) = rec(15).data;
      params.pass.SiderealTime(:,I) = rec(16).data;
      params.pass.nSPI(:,I) = rec(17).data;
      params.pass.resampleRate1(:,I) = rec(18).data;
      params.pass.resampleRate2(:,I) = rec(19).data;
      params.pass.nResample1(:,I) = rec(20).data;
      params.pass.nResample2(:,I) = rec(21).data;
      params.pass.betaParam1(:,I) = rec(22).data;
      params.pass.betaParam2(:,I) = rec(23).data;
    end % doStochmap if statement

  end % badtimesflag if statement

end % loop over intervals I

if doStochmap
  params.frompreproc = true;
  params.jobNumber = jobNumber;
  params = addSTAMPparams(params);
  fprintf('Done with preproc phase, beginning stochmap phase...\n');

  % check times and calculate beg_pause and end_pause
  [startGPS, endGPS, params] = check_times(params);

  if bknd_study==1      
   available_segments=params.pass.GPStimes;

    fid=fopen(['results/' num2str(startGPS) '.txt'],'w+');
    BkndTrials = floor((endGPS-startGPS)/bknd_study_dur);
    for ii=1:BkndTrials
      fprintf('trial = %i\n',ii);

      % bstart is the gps time for the first segment for this background trial
      bstart = startGPS + (ii-1)*bknd_study_dur;
      % bend is the gps time for last segment for this background trial
      bend = bstart + bknd_study_dur;
      
      params.pass.GPStimes = bstart:params.pass.segmentDuration/2:bend;

      % which_beg_idx is the gps time that is closest to the requested gps 
      % time for this background trial.  (This is necessary to handle the 
      % situation where the requested time is not a segment start time.)
      which_beg_idx=min(abs(available_segments-bstart));
      beg_idx = min(find(abs(available_segments-bstart)==which_beg_idx));
      which_end_idx=min(abs(available_segments-bend));
      end_idx = max(find(abs(available_segments-bend)==which_end_idx));
      params.pass.which_segs=beg_idx:1:end_idx;

      stoch_out = stochmap(params,jobsFile,bstart,bend);
      if ii==1
        fprintf(fid, '%% GPS ');
        for ss=1:length(stoch_out.search)
          fprintf(fid,'%s ', stoch_out.search{ss});
        end
        fprintf(fid,'\n');
      end
   
      fprintf(fid,'%f ',bstart);
      for ss=1:length(stoch_out.max_SNR)
        fprintf(fid, '%f ', stoch_out.max_SNR(ss));
      end
      fprintf(fid, '\n');
    end
    fclose(fid);
    preproc_out=stoch_out;
  else
    stochmap(params,jobsFile,startGPS,endGPS);
    preproc_out=0;
  end
end

% Define preproc_out.done so that the output argument is always defined.
preproc_out.done = true;
fprintf('elapsed_time = %f\n', toc(tStart));

return
