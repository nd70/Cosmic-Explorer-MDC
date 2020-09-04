% STRAINCHECK(paramsFile, jobsFile)
%
% The aim if this script is to compare the band-limited RMS of two channels.
% Usually the two channels would be AS_Q (after applying the response function)
% and the calibrated strain. The script goes through the list of science
% segments given in jobsFile and calculates the band-limited RMS of the
% channels specified in the paramsFile, and the RMS of their difference.
% These values are stored as the time-series rms1, rms2 and rmsDiff in a
% location specified by the paramsFile.
%
function straincheck(paramsFile, jobsFile)

% read in params structure from a file
params = readParamsFromFile(paramsFile);

% extract parameters from the params structure setting to default values
% if necessary
try
  suppressFrWarnings = params.suppressFrWarnings;
  if (suppressFrWarnings)
    warning off frgetvect:info;
  else
    warning on frgetvect:info;
  end;
catch
  % Do nothing
end;
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
  segmentDuration = params.segmentDuration;
catch
  error('segmentDuration parameter not set');
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
  resampleRate = params.resampleRate1;
catch
  error('resampleRate parameter not set');
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
  error('frameDuration1 parameter not set')
end
try
  frameDuration2 = params.frameDuration2;
catch
  error('frameDuration2 parameter not set')
end
try
  nResample1 = params.nResample1;
catch
  nResample1 = 10;
end
try
  nResample2 = params.nResample2;
catch
  nResample2 = 10;
end
try
  betaParam1 = params.betaParam1;
catch
  betaParam1 = 5;
end
try
  betaParam2 = params.betaParam2;
catch
  betaParam2 = 5;
end
try
  highPassFreq1 = params.highPassFreq1;
catch
  highPassFreq1 = 40;
end
try
  highPassFreq2 = params.highPassFreq2;
catch
  highPassFreq2 = 40;
end
try
  highPassOrder1 = params.highPassOrder1;
catch
  highPassOrder1 = 6;
end
try
  highPassOrder2 = params.highPassOrder2;
catch
  highPassOrder2 = 6;
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
  outputFilePrefix = params.outputFilePrefix;
catch
  error('outputFilePrefix parameter not set');
end
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

% Padding for filters and downsampling
nPad = 16384;

kappa = 1.0E+020;

% set values for data windowing, zero-padding, and FFT
numPoints = segmentDuration*resampleRate; 
dataWindow = hanning(numPoints);
%dataWindow = tukeywin(numPoints, 0.5);
%dataWindow = ones(numPoints, 1);
dataWindow = sqrt(numPoints)*dataWindow/norm(dataWindow);
fftLength = numPoints;

% Frequency spacing
deltaF = resampleRate/fftLength;

% read in start times and durations from the job file
[startTimes, jobDurations] = readJobsFile(jobsFile, params.jobsFileCommentStyle);

% Make a sorted list of science segments
segments = [ startTimes, jobDurations ];
segments = sortrows(segments);
scienceRunStartTime = segments(1, 1);
scienceRunEndTime = segments(end, 1) + segments(end, 2);
totalDuration = scienceRunEndTime - scienceRunStartTime;

% Set up the frequency indices
flowIdx = floor(flow/deltaF);
flow = flowIdx*deltaF;
fhighIdx = ceil(fhigh/deltaF);
fhigh = fhighIdx*deltaF;
numFreqs = fhighIdx - flowIdx + 1;
fIdx = [flowIdx+1:fhighIdx+1];

% construct filter coefficients for high-pass filtering
if doHighPass1
  [b1, a1] = butter(highPassOrder1, highPassFreq1/(resampleRate/2), 'high');
end
if doHighPass2
  [b2, a2] = butter(highPassOrder2, highPassFreq2/(resampleRate/2), 'high');
end

% construct frequency mask for later use
mask.data = constructFreqMask(flow, fhigh, deltaF, ...
                              freqsToRemove, nBinsToRemove, doFreqMask);
mask.flow = flow;
mask.deltaF = deltaF;

% Chi-squared degree of the sum of |h_k|^2 is given by the 2 x number
% of independent bins
nBins = sum(mask.data);

% channel names
channelName1 = [ ifo1 ':' ASQchannel1 ];
channelName2 = [ ifo2 ':' ASQchannel2 ];

% common max frame duration
frameDuration = max([ frameDuration1, frameDuration2 ]);

% How many segments per frame? must divide evenly
if (mod(frameDuration, segmentDuration) ~= 0)
  error('Segments must evenly divide frame length');
end;
segsPerFrame = frameDuration/segmentDuration;

% read in calibration info
if ( ~strncmp(alphaBetaFile1,   'none', length(alphaBetaFile1))   & ...
     ~strncmp(calCavGainFile1,  'none', length(calCavGainFile1))  & ...
     ~strncmp(calResponseFile1, 'none', length(calResponseFile1)) )
  useCal1 = true;
  [t1, f1, R01, C01, alpha1, gamma1] = ...
    readCalibrationFromFiles(alphaBetaFile1, calCavGainFile1, calResponseFile1);
else
  useCal1 = false;
end;

if ( ~strncmp(alphaBetaFile2,   'none', length(alphaBetaFile2))   & ...
     ~strncmp(calCavGainFile2,  'none', length(calCavGainFile2))  & ...
     ~strncmp(calResponseFile2, 'none', length(calResponseFile2)) )
  useCal2 = true;
  [t2, f2, R02, C02, alpha2, gamma2] = ...
    readCalibrationFromFiles(alphaBetaFile2, calCavGainFile2, calResponseFile2);
else
  useCal2 = false;
end;

% Loop over the list of science segments
numSegments = size(segments, 1);
for segment = 1:numSegments

  jobDuration = segments(segment, 2);

  % Span of the current science segment
  startTime = segments(segment, 1);
  endTime = startTime + jobDuration;

  fprintf('Science segment %d/%d: %d-%d', segment, numSegments, ...
    startTime, startTime + jobDuration);

  % Ignore science segments less than frameDuration long
  if (jobDuration < frameDuration)
    fprintf(' (too short, skipping)\n');
    continue;
  end;
  fprintf('\n');
  
  % Get the cache info for this segment, skipping cases where it is not available
  [gpsList1, fileList1, status, durList1] = ...
    framelist(ifo1(1), frameType1, startTime, jobDuration);
  if (status ~= 0)
    fprintf('straincheck: error getting frame-file list for frame data: %s, %s, %d-%d\n', ifo1(1), frameType1, startTime, startTime+jobDuration);
    continue;
  end;
  
  [gpsList2, fileList2, status, durList2] = ...
    framelist(ifo2(1), frameType2, startTime, jobDuration);
  if (status ~= 0)
    fprintf('straincheck: error getting frame-file list for frame data: %s, %s, %d-%d\n', ifo2(1), frameType2, startTime, startTime+jobDuration);
    continue;
  end;
  
  % Create the time-series that will store the results of this analysis
  % Note that in the current implementation we will miss out on any partial
  % frames at the end
  nFrames = floor(jobDuration/frameDuration);
  nOutput = nFrames*segsPerFrame;
  result{segment}.rms1 = constructTimeSeries(repmat(NaN, 1, nOutput), startTime, ...
    segmentDuration, 0.0, 0.0);
  result{segment}.rms2 = result{segment}.rms1;
  result{segment}.rmsDiff = result{segment}.rms1;

  % Loop over the subframes of the current science segment. Note that the data is
  % read in chunks the size of a frame (whatever that is set to be) but probably
  % won't line up with frame-file boundaries
  rmsIdx = 0;
  intervalStartTime = startTime;
  while (intervalStartTime + frameDuration <= endTime)

    fprintf('  Frame interval %d-%d\n', ...
      intervalStartTime, intervalStartTime + frameDuration);

    [adcdata1, dataOK] = readData(channelName1, intervalStartTime, ...
      frameDuration, frameType1, gpsList1, fileList1, durList1);
  
    if (~dataOK)
      fprintf('Bad data for %s at %d\n', ifo1, intervalStartTime);
      rmsIdx = rmsIdx + segsPerFrame;
      intervalStartTime = startTime + rmsIdx*segmentDuration;
      continue;
    end;
  
    [adcdata2, dataOK] = readData(channelName2, intervalStartTime, ...
      frameDuration, frameType2, gpsList2, fileList2, durList2);

    if (~dataOK)
      fprintf('Bad data for %s at %d\n', ifo2, intervalStartTime);
      rmsIdx = rmsIdx + segsPerFrame;
      intervalStartTime = startTime + rmsIdx*segmentDuration;
      continue;
    end;
  
    %%
    %% Pre-process IFO 1 data
    %%

    % Pad with nPad samples at each end
    adcdata1.data = [ zeros(nPad, 1); adcdata1.data; zeros(nPad, 1) ];

    % Downsample
    sampleRate1 = 1/adcdata1.deltaT;
    if (sampleRate1 ~= resampleRate)
      adcdata1.data = resample(adcdata1.data, 1, ...
        floor(sampleRate1/resampleRate), nResample1, betaParam1);
      adcdata1.deltaT = 1/resampleRate;
    end
  
    % High-pass filter
    if doHighPass1
      adcdata1.data = filtfilt(b1, a1, adcdata1.data);
    end

    % Unpad
    nPad2 = floor(nPad*resampleRate/sampleRate1);
    adcdata1.data = adcdata1.data(nPad2:end-nPad2-1);

    %%
    %% Pre-process IFO 2 data
    %%

    % Scale strain data
    adcdata2.data = kappa*adcdata2.data;

    % Pad with nPad samples at each end
    adcdata2.data = [ zeros(nPad, 1); adcdata2.data; zeros(nPad, 1) ];

    % Downsample
    sampleRate2 = 1/adcdata2.deltaT;
    if (sampleRate2 ~= resampleRate)
      adcdata2.data = resample(adcdata2.data, 1, ...
        floor(sampleRate2/resampleRate), nResample2, betaParam2);
      adcdata2.deltaT = 1/resampleRate;
    end;
  
    % High-pass
    if doHighPass2
      adcdata2.data = filtfilt(b2, a2, adcdata2.data);
    end
  
    % Unpad
    nPad2 = floor(nPad*resampleRate/sampleRate2);
    adcdata2.data = adcdata2.data(nPad2:end-nPad2-1);

    % Calculate spectra
    [P1, fs1] = psd(adcdata1.data, fftLength, resampleRate, [], fftLength/2);
    P1 = P1(fIdx);
    fs1 = fs1(fIdx);
    P2 = psd(adcdata2.data, fftLength, resampleRate, [], fftLength/2);
    P2 = P2(fIdx);

    % For each subinterval within a frame
    for m = 0:segsPerFrame-1

      fprintf('    Subinterval %d/%d:%d-%d\n', m+1, segsPerFrame, ...
        intervalStartTime, intervalStartTime + segmentDuration);

      % calculate response functions from calibration data
      if (useCal1)
        [R1, responseOK] = calculateResponse(t1, f1, R01, C01, alpha1, ...
                             gamma1, intervalStartTime + segmentDuration/2, ASQchannel1);
        if (~responseOK)
          fprintf('Bad response function for %s at %d\n', ifo1, intervalStartTime);
          rmsIdx = rmsIdx + 1;
          intervalStartTime = startTime + rmsIdx*segmentDuration;
          continue;
        end
  
        % evaluate response function at desired frequencies
        resp1 = convertResponse(f1, kappa*R1, flow, deltaF, numFreqs, 0, 0);
      else
        % the data is already calibrated
        resp1 = constructFreqSeries(kappa*ones(numFreqs, 1), flow, deltaF);
      end
  
      if (useCal2)
        [R2, responseOK] = calculateResponse(t2, f2, R02, C02, alpha2, ...
                             gamma2, intervalStartTime + segmentDuration/2, ASQchannel2);
        if (~responseOK)
          fprintf('Bad response function for %s at %d\n', ifo2, intervalStartTime);
          rmsIdx = rmsIdx + 1;
          intervalStartTime = startTime + rmsIdx*segmentDuration;
          continue;
        end
  
        % evaluate response function at desired frequencies
        resp2 = convertResponse(f2, R2, flow, deltaF, numFreqs, 0, 0);
      else
        % the data is already calibrated
        resp2 = constructFreqSeries(ones(numFreqs, 1), flow, deltaF);
      end

      % Work out the new spectra after convolving with the response functions
      % NOTE: probably not a good idea for the spectrum from DARM_ERR, since we
      % would be convolving the DARM_ERR spectrum with a "local" version of R(f)
      % specific to this time interval. Instead I will only use the spectrum
      % as calculated from calibrated strain and find differences with respect
      % to that.
      PP1 = abs(resp1.data).^2.*P1;
      PP2 = abs(resp2.data).^2.*P2;

      % Window, fft, apply response function and mask
      rtilde1 = fft(adcdata1.data(m*numPoints+1:(m+1)*numPoints).*dataWindow, fftLength);
      htilde1 = mask.data.*resp1.data.*rtilde1(fIdx);

      rtilde2 = fft(adcdata2.data(m*numPoints+1:(m+1)*numPoints).*dataWindow, fftLength);
      htilde2 = mask.data.*resp2.data.*rtilde2(fIdx);

      % |h1[k]|^2 - the weighting is meant to ensure that the series is
      % distributed as a chi-square of degree 2
      hsq1 = (2/fftLength)*abs(htilde1).^2./PP2;

      % |h2[k]|^2 wrt to spectrum PP2
      hsq2 = (2/fftLength)*abs(htilde2).^2./PP2;

      % |h2[k] - h1[k]|^2 wrt to spectrum PP2  
      hdiffsq = (2/fftLength)*abs(htilde2 - htilde1).^2./PP2;

      % Since each htilde is weighted by the inverse of the spectrum, with the additional factor
      % of 2/N, |h_k|^2 should be distributed like chi^2 degree 2. Assuming each
      % frequency is independent, the sum over the frequency range is thus chi^2 degree
      % 2*sum(mask.data) (since we don't count bins that are is masked out). This degree
      % will be large, so the distribution of each rms will be approximately normal with
      % mean 2*sum(mask.data), variance 4*sum(mask.data)
      r1 = sum(hsq1);
      r2 = sum(hsq2);
      rd = sum(hdiffsq);

      result{segment}.rms1.data(rmsIdx+1) = r1;
      result{segment}.rms2.data(rmsIdx+1) = r2;
      result{segment}.rmsDiff.data(rmsIdx+1) = rd;

      % For debugging
      z = 100*sqrt(rd/r1);

      %plot2resp(resp1);
      %plot2time(adcdata1, adcdata2);
      %plot3freq(htilde1, htilde2);

      rmsIdx = rmsIdx + 1;
      intervalStartTime = startTime + rmsIdx*segmentDuration;

    end; % loop over subintervals
  end; % loop over frames

  save([outputFilePrefix 'straincheck.mat'], 'result', 'nBins');

end; % loop over science segments

return; % End of main function straincheck

function plot2resp(resp)

  f = resp.flow + resp.deltaF*[0:length(resp.data)-1];
  subplot(2,1,1);
  plot(f, abs(resp.data));
  subplot(2,1,2);
  plot(f, angle(resp.data));
  pause;

return;

function plot2time(adcdata1, adcdata2)

  subplot(2, 1, 1);
  plot(adcdata1.data(100:end-50));
  subplot(2, 1, 2);
  plot(adcdata2.data(100:end-50));
  pause;
  
return;

function plot3freq(htilde1, htilde2)

  f = htilde1.flow + htilde1.deltaF*[0:length(htilde1.data)-1];
  %plot(abs(resp2.data));
  %plot(f, (angle(resp1.data)));
  subplot(3, 1, 1);
  plot(f, (abs(htilde1.data)));
  subplot(3, 1, 2);
  plot(f, (abs(htilde2.data)));
  subplot(3, 1, 3);
  plot(f, abs(htilde1.data) - abs(htilde2.data));
  %plot(f, 2*(abs(htilde2.data) - abs(htilde1.data))/(rms2.data(rmsIdx) + ...
  %  rms1.data(rmsIdx)));
  %plot(f, (angle(htilde1.data))-(angle(htilde2.data)), 'r');
  pause;
  
return;

function [adcdata, dataOK] = readData(channelName, dataStartTime, ...
                       dataDuration, frameType, gpsList, fileList, durList)

  % initially set dataOK to true
  dataOK = true;

  % get the data
  chanObject  = chanstruct(channelName);
  chanObject.type = frameType;
  [vector, sampleRate, vectorError] ...
    = chanvector(chanObject, dataStartTime, dataDuration, gpsList, ...
        fileList, durList);

  % check that the data is OK
  if (vectorError == 0)
    dataOK = true;
  else
    fprintf('readData: missing frame data: %s, %s, %d-%d\n', ...
      channelName, frameType, dataStartTime, dataStartTime+dataDuration);
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
  end;

return; % readData
