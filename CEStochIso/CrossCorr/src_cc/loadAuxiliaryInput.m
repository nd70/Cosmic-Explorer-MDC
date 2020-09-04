function params=loadAuxiliaryInput(params)

% loads and initializes various auxiliary data
% Previously done in the old stochastic.m
% 
% input and output: params struct
%
% Routine copied from stochastic.m and modified my Stefan Ballmer 
% sballmer@caltech.edu 
%
% $Id: loadAuxiliaryInput.m,v 1.19 2009-02-10 05:11:21 ethrane Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get appropriate detector structure for each ifo
if (isnan(params.azimuth1))
  params.detector1 = getdetector(params.site1);
else
  params.detector1 = getdetector(params.site1,params.azimuth1);
end
if (isnan(params.azimuth2))
  params.detector2 = getdetector(params.site2);
else
  params.detector2 = getdetector(params.site2,params.azimuth2);
end

% construct filter coefficients for high-pass filtering
if params.doHighPass1
  [params.filt1.b,params.filt1.a] = butter(params.highPassOrder1, params.highPassFreq1/(params.resampleRate1/2), 'high');
end

if params.doHighPass2
  [params.filt2.b,params.filt2.a] = butter(params.highPassOrder2, params.highPassFreq2/(params.resampleRate2/2), 'high');
end

% set values for psd estimation (on resampled data, HP filtered data)
params.psd1.FFTLength = params.resampleRate1*(1/params.deltaF);
params.psd1.Window    = hann(params.psd1.FFTLength);
params.psd1.OverlapLength = params.psd1.FFTLength/2; 
params.psd1.detrendFlag  = 'none';

params.psd2.FFTLength = params.resampleRate2*(1/params.deltaF);
params.psd2.Window    = hann(params.psd2.FFTLength);
params.psd2.OverlapLength = params.psd2.FFTLength/2; 
params.psd2.detrendFlag  = 'none';

% set values for data windowing, zero-padding, and FFT
if params.doOverlap
  params.hannDuration1=params.segmentDuration;
  params.hannDuration2=params.segmentDuration;
end;
numPoints1    = params.segmentDuration*params.resampleRate1; 
params.fft1.dataWindow   = tukeywin(numPoints1, params.hannDuration1/params.segmentDuration);
params.fft1.fftLength    = 2*numPoints1;

numPoints2    = params.segmentDuration*params.resampleRate2;
params.fft2.dataWindow   = tukeywin(numPoints2, params.hannDuration2/params.segmentDuration);
params.fft2.fftLength    = 2*numPoints2;

% construct frequency mask for later use
data = constructFreqMask(params.flow, params.fhigh, params.deltaF, ...
                         params.freqsToRemove, params.nBinsToRemove, params.doFreqMask);
params.mask = constructFreqSeries(data, params.flow, params.deltaF);

% calculate overlap reduction function
params.f = params.flow + params.deltaF*transpose([0:params.numFreqs-1]);
if params.doDirectional
  % gamma equal one - correction for gamma0 and tau is done in ccStatReadout
  params.gamma = constructFreqSeries(ones(params.numFreqs,1), params.flow, params.deltaF, 1);
  if params.doAllSkyComparison
    % also calculate the All-Sky overlap reduction function
    data  = overlapreductionfunction(params.f, params.detector1, params.detector2);
    % If one of the data streams is params.heterodyned, set params.gamma.symmetry to 0
    % so that negative frequencies are not included in the normalization.
    params.gammaAllSky = constructFreqSeries(data, params.flow, params.deltaF, 1);
  end
else
  data  = overlapreductionfunction(params.f, params.detector1, params.detector2);
  % If one of the data streams is params.heterodyned, set params.gamma.symmetry to 0
  % so that negative frequencies are not included in the normalization.
  if params.heterodyned
    params.gamma = constructFreqSeries(data, params.flow, params.deltaF, 0);
  else
    params.gamma = constructFreqSeries(data, params.flow, params.deltaF, 1);
  end
end

% construct name of gps times files and frame cache files for this job
if ~params.intermediate %cet-----------------------------------------------------
  if (strcmp(params.gpsTimesPath1, '<auto>'))
    params.gpsTimesFile1 = '<auto>';
  else
    params.gpsTimesFile1 = ...
        [params.gpsTimesPath1 'gpsTimes' params.ifo1(1) '.' num2str(params.jobNumber) '.txt'];
  end;

  if (strcmp(params.gpsTimesPath2, '<auto>'))
    params.gpsTimesFile2 = '<auto>';
  else
    params.gpsTimesFile2 = ...
        [params.gpsTimesPath2 'gpsTimes' params.ifo2(1) '.' num2str(params.jobNumber) '.txt'];
  end;

  params.frameCacheFile1 = ...
    [params.frameCachePath1 'frameFiles' params.ifo1(1) '.' num2str(params.jobNumber) '.txt'];
  params.frameCacheFile2 = ...
    [params.frameCachePath2 'frameFiles' params.ifo2(1) '.' num2str(params.jobNumber) '.txt'];
else
%cetfeb9  params.intFrameFile = ...
%cetfeb9   [params.intFrameCachePath 'frameFiles' params.ifo1(1) '.' num2str(params.jobNumber) '.txt'];
  params.intFrameFile = ...
    [params.intFrameCachePath 'frameFiles' params.ifo1(1) params.ifo2(1) '.' num2str(params.jobNumber) '.txt'];
end %--------------------------------------------------------------------------
if params.doInjFromFile1
  if (strcmp(params.injGPSTimesPath1, '<auto>'))
    params.injGPSTimesFile1 = '<auto>';
  else
    params.injGPSTimesFile1 = ...
        [params.injGPSTimesPath1 'gpsTimes' params.ifo1(1) '.' num2str(params.jobNumber) '.txt'];
  end;
  params.injFrameCacheFile1 = ...
      [params.injFrameCachePath1 'frameFiles' params.ifo1(1) '.' num2str(params.jobNumber) '.txt'];
end
if params.doInjFromFile2
  if (strcmp(params.injGPSTimesPath2, '<auto>'))
    params.injGPSTimesFile2 = '<auto>';
  else
    params.injGPSTimesFile2 = ...
        [params.injGPSTimesPath2 'gpsTimes' params.ifo2(1) '.' num2str(params.jobNumber) '.txt'];
  end;
  params.injFrameCacheFile2 = ...
      [params.injFrameCachePath2 'frameFiles' params.ifo2(1) '.' num2str(params.jobNumber) '.txt'];
end

% channel names
params.channelName1 = [params.ifo1 ':' params.ASQchannel1];
params.channelName2 = [params.ifo2 ':' params.ASQchannel2];
if params.doInjFromFile1
  params.injChannelName1 = [params.injPrefix1 ':' params.injChannel1];
end
if params.doInjFromFile2
  params.injChannelName2 = [params.injPrefix2 ':' params.injChannel2];
end

% read in calibration info

if ( ~strncmp(params.alphaBetaFile1,   'none', length(params.alphaBetaFile1))   & ...
     ~strncmp(params.calCavGainFile1,  'none', length(params.calCavGainFile1))  & ...
     ~strncmp(params.calResponseFile1, 'none', length(params.calResponseFile1)) )
[params.cal1.t, params.cal1.f, params.cal1.R0, params.cal1.C0, params.cal1.alpha, params.cal1.gamma] = ...
  readCalibrationFromFiles(params.alphaBetaFile1, params.calCavGainFile1, params.calResponseFile1);
  params.channel1Calibrated = false;
else
  params.cal1 = getTrivialCalibration();
  params.channel1Calibrated = true;
end;

if ( ~strncmp(params.alphaBetaFile2,   'none', length(params.alphaBetaFile2))   & ...
     ~strncmp(params.calCavGainFile2,  'none', length(params.calCavGainFile2))  & ...
     ~strncmp(params.calResponseFile2, 'none', length(params.calResponseFile2)) )
[params.cal2.t, params.cal2.f, params.cal2.R0, params.cal2.C0, params.cal2.alpha, params.cal2.gamma] = ...
  readCalibrationFromFiles(params.alphaBetaFile2, params.calCavGainFile2, params.calResponseFile2);
  params.channel2Calibrated = false;
else
  params.cal2 = getTrivialCalibration();
  params.channel2Calibrated = true;
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initiate params.SkyPattern
if params.doDirectional
  if params.useSkyPatternFile
  params.SkyPattern=load(params.SkyPatternFile); % SkyPatternFile should be in hours and degrees for RA & declination, respectively
  else
    params.SkyPattern=SkyPattern(params.SkyPatternRightAscensionNumPoints,params.SkyPatternDeclinationNumPoints);
  end
end
 
% load signal power specrum file Hf
if params.useSignalSpectrumHfFromFile
  Hfdata=load(params.HfFile);
  if params.HfFileInterpolateLogarithmic
    params.Hf=constructFreqSeries(exp(interp1(log(Hfdata(:,1)),log(Hfdata(:,2)),log(params.flow:params.deltaF:params.fhigh)))', params.flow, params.deltaF, 1);
  else
    params.Hf=constructFreqSeries(interp1(Hfdata(:,1),Hfdata(:,2),params.flow:params.deltaF:params.fhigh)', params.flow, params.deltaF, 1);
  end
  clear Hfdata;
end
% if in radiometer mode params.fRef contains the power spectrum Hf
if params.doDirectional
  params.fRef=params.Hf;
  params.alphaExp=[];
end


% load the simulated point source related files
if params.doSimulatedPointSource
  filename=[params.simulationPath,'/',params.simulatedPointSourcesFile];
  simPSSkyPositions=load(filename);
  filename=[params.simulationPath,'/',params.simulatedPointSourcesPowerSpec];
  simPSPowerSpec=load(filename);
  % prepare the point source simulation
  ra=simPSSkyPositions(:,1);
  decl=simPSSkyPositions(:,2);
  if size(simPSSkyPositions,2)>2
    power=simPSSkyPositions(:,3);
  else
    power=ones(size(ra));
  end
  initPointSourceData(params.resampleRate1,simPSPowerSpec,params.simulatedPointSourcesInterpolateLogarithmic,...
                      params.detector1,params.detector2,ra,decl,power,params.flow,params.deltaF,params.numFreqs,...
                      params.cal1.t,params.cal1.f,params.cal1.R0,params.cal1.C0,params.cal1.alpha,params.cal1.gamma,...
		      params.ASQchannel1,params.alphaBetaFile1,params.calCavGainFile1,params.calResponseFile1,...
		      params.cal2.t,params.cal2.f,params.cal2.R0,params.cal2.C0,params.cal2.alpha,params.cal2.gamma,...
		      params.ASQchannel2,params.alphaBetaFile2,params.calCavGainFile2,params.calResponseFile2,...
		      params.simulatedPointSourcesBufferDepth,params.simulatedPointSourcesHalfRefillLength,...
		      params.simulatedPointSourcesMakeIncoherent);
end;

% load the simulated sky map related files
if params.doSimulatedSkyMap
  filename=[params.simulationPath,'/',params.simulatedSkyMapFile];
  [map,coord1,coord2]=readMap(filename,params.simulatedSkyMapFileType,...
                                       params.simulatedSkyMapFileNumBins,...
				       params.simulatedSkyMapInjectAsSpH,...
				       params.simulatedSkyMapConvertLmax,...
				       params.simulatedSkyMapConvertDeg);
  filename=[params.simulationPath,'/',params.simulatedSkyMapPowerSpec];
  simSkyPowerSpec=load(filename);
  initSkyMapData(params.simulatedSkyMapInjectTimeDomain,...
                 params.resampleRate1,params.resampleRate2,...
                 params.nResample1,params.nResample2,params.betaParam1,params.betaParam2,...
                 simSkyPowerSpec,params.simulatedSkyMapInterpolateLogarithmic,...
                 params.simulatedSkyMapConvertLmax, params.ifo1, params.ifo2, params.gammaLM_coeffsPath,...
                 params.detector1,params.detector2,params.simulatedSkyMapInjectAsSpH,coord1,coord2,map,...
		 params.flow,params.deltaF,params.numFreqs,...
                 params.cal1.t,params.cal1.f,params.cal1.R0,params.cal1.C0,params.cal1.alpha,params.cal1.gamma,...
		 params.ASQchannel1,params.alphaBetaFile1,params.calCavGainFile1,params.calResponseFile1,...
		 params.cal2.t,params.cal2.f,params.cal2.R0,params.cal2.C0,params.cal2.alpha,params.cal2.gamma,...
		 params.ASQchannel2,params.alphaBetaFile2,params.calCavGainFile2,params.calResponseFile2,...
		 params.simulatedSkyMapBufferDepth,params.simulatedSkyMapHalfRefillLength,...
		 params.simulatedSkyMapMakeIncoherent);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the simulated detector noise files
if params.doSimulatedDetectorNoise
  filename=[params.simulationPath,'/',params.simulatedDetectorNoisePowerSpec1];
  simDetectorNoisePowerSpec1=load(filename);
  filename=[params.simulationPath,'/',params.simulatedDetectorNoisePowerSpec2];
  simDetectorNoisePowerSpec2=load(filename);
  initDetectorNoiseData(params.resampleRate1,params.resampleRate2,...
                 params.nResample1,params.nResample2,params.betaParam1,params.betaParam2,...
                 simDetectorNoisePowerSpec1,simDetectorNoisePowerSpec2,...
                 params.simulatedDetectorNoiseInterpolateLogarithmic,...
		 params.flow,params.deltaF,params.numFreqs,...
                 params.cal1.t,params.cal1.f,params.cal1.R0,params.cal1.C0,params.cal1.alpha,params.cal1.gamma,...
		 params.ASQchannel1,params.alphaBetaFile1,params.calCavGainFile1,params.calResponseFile1,...
		 params.cal2.t,params.cal2.f,params.cal2.R0,params.cal2.C0,params.cal2.alpha,params.cal2.gamma,...
		 params.ASQchannel2,params.alphaBetaFile2,params.calCavGainFile2,params.calResponseFile2,...
		 params.simulatedDetectorNoiseBufferDepth,params.simulatedDetectorNoiseHalfRefillLength);
end;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set number of trials to 1 if no monte carlo simulations
if params.doMonteCarlo==false; 
%  August 3, 2010: Overriding the numTrials=1 default so that I can try using
%  multiple trials to do other things, such as looping over time shifts.
  params.numTrials = 1; 
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load timing correction
if params.doTimingTransientSubtraction1
  params.TimingTransient1=load(params.TimingTransientFile1);
end
if params.doTimingTransientSubtraction2
  params.TimingTransient2=load(params.TimingTransientFile2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check that params.numSegmentsPerInterval is odd    
if mod(params.numSegmentsPerInterval,2)==0
  error('params.numSegmentsPerInterval must be odd');
end

% determine number of intervals and segments to analyse
params.bufferSecsMax = max(params.bufferSecs1,params.bufferSecs2);

params.M = floor( (params.jobDuration - 2*params.bufferSecsMax)/params.segmentDuration );

if ~params.doSidereal %cet---------------------------------------------------
  if params.doOverlap
    params.numSegmentsTotal = 2*params.M-1;
    params.numIntervalsTotal = 2*(params.M - ...
				  (params.numSegmentsPerInterval-1)) - 1;
    params.intervalTimeStride = params.segmentDuration/2;
  else 
    params.numSegmentsTotal = params.M;
    params.numIntervalsTotal = params.M - (params.numSegmentsPerInterval-1);
    params.intervalTimeStride = params.segmentDuration;
  end
  params.centeredStartTime = params.startTime + params.bufferSecsMax + ...
    floor( (params.jobDuration - 2*params.bufferSecsMax - ...
	    params.M*params.segmentDuration)/ 2 );
else %calculate parameters compatible with the sidereal time------------------
  srfac = 23.9344696 / 24; %sidereal time conversion factor
  startTime = params.startTime;
  srtime = GPStoGreenwichMeanSiderealTime(startTime) * 3600;
  md = mod(srtime,params.segmentDuration/srfac);
  params.centeredStartTime = round(startTime + ...
				   params.segmentDuration - md*srfac);
  if params.centeredStartTime - startTime < params.bufferSecsMax
    params.centeredStartTime = params.centeredStartTime + ...
    params.segmentDuration;
  end
  %now we can calculate the remaining bookkeeping variables
  params.M = floor( (params.jobDuration - params.bufferSecsMax - ...
		     params.centeredStartTime + ...
		     params.startTime) / params.segmentDuration );
  if params.doOverlap
    params.numIntervalsTotal = 2*(params.M - ...
				  (params.numSegmentsPerInterval-1)) - 1;
    params.intervalTimeStride = params.segmentDuration/2;
  else
    params.numIntervalsTotal = params.M - (params.numSegmentsPerInterval-1);
    params.intervalTimeStride = params.segmentDuration;
  end
end %--------------------------------------------------------------------------

if params.intermediate %cet----------------------------------------------------
  params.intFrameFile;
  intFrameList = textread(params.intFrameFile,'%s');
  [listsize,itemsize] = size(intFrameList);
  params.numIntervalsTotal = listsize;
  params.intervalTimeStride = params.segmentDuration/2;
end %--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.lastLoadedDataEnd1 = 0;
params.lastLoadedDataEnd2 = 0;

