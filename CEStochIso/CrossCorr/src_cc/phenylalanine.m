function [calPSD, rbartilde] = phenylalanine(params, gps, data_dur, ...
					     fr_gps, fr_url, fr_durs, ifo)
% by Eric Thrane on August 25, 2010
% phenylalanine is an ingredient in stochastic_lite.
% It calculates calPSD and rbartilde for a given segment and for a given ifo.

  % select parameters based on value of ifo
  if ifo==1
    params.bufferSecs = params.bufferSecs1;
    params.channelName = params.channelName1;
    params.resampleRate = params.resampleRate1;
    params.betaParam = params.betaParam1;
    params.doHighPass = params.doHighPass1;
    params.nResample = params.nResample1;
    params.filt.a = params.filt1.a;
    params.filt.b = params.filt1.b;
    params.PSD.FFTLength = params.psd1.FFTLength;
    params.PSD.Window =  params.psd1.Window;
    params.PSD.OverlapLength =  params.psd1.OverlapLength;
    params.PSD.detrendFlag =  params.psd1.detrendFlag;
    params.fft.dataWindow = params.fft1.dataWindow;
    params.fft.fftLength = params.fft1.fftLength;
  else
    params.bufferSecs = params.bufferSecs2;
    params.channelName = params.channelName2;
    params.resampleRate = params.resampleRate2;
    params.betaParam = params.betaParam2;
    params.doHighPass = params.doHighPass2;
    params.nResample = params.nResample2;
    params.filt.a = params.filt2.a;
    params.filt.b = params.filt2.b;
    params.PSD.FFTLength = params.psd2.FFTLength;
    params.PSD.Window =  params.psd2.Window;
    params.PSD.OverlapLength =  params.psd2.OverlapLength;
    params.PSD.detrendFlag =  params.psd2.detrendFlag;
    params.fft.dataWindow = params.fft2.dataWindow;
    params.fft.fftLength = params.fft2.fftLength;
  end

  % construct chanObjects for call to chanvector
  chanObject  = chanstruct(params.channelName);

  % read in time series for each segment
  [vector, sampleRate, vectorError] = chanvector(chanObject, ...
    gps, data_dur, fr_gps, fr_url, fr_durs);

  % check that the data is OK
  if vectorError == 0
    % fill time-series data structures
    adcdata.data   = transpose(vector);
    adcdata.tlow   = gps;
    adcdata.deltaT = 1/sampleRate;
  else
    fprintf('Error reading in data.\n');
  end

  % create fbase and phase fields
  adcdata.fbase = NaN;
  adcdata.phase = NaN;

  % downsample the data
  sampleRate = 1/adcdata.deltaT;
  p = 1;
  q = floor(sampleRate/params.resampleRate);
  deltaT = 1/params.resampleRate;

  if sampleRate == params.resampleRate
    data = adcdata.data;
  else
    data = resample(adcdata.data, p, q, params.nResample, params.betaParam);
  end

  % if requested overwrite data with simulated detector noise
  % simulated detector noise is already resampled
  if params.doSimulatedDetectorNoise
    % inject white noise with RMS=1
    data = randn(size(data));
  end

  o = constructTimeSeries(data, adcdata.tlow, deltaT, adcdata.fbase, ...
			  adcdata.phase);

  % free-up some memory
  clear adcdata

  % high-pass the data
  if params.doHighPass
    highpassed = constructTimeSeries(filtfilt(params.filt.b, ...
      params.filt.a, o.data), o.tlow, o.deltaT, o.fbase, o.phase);
  else
    highpassed = o;
  end

  % chop-off bad data at start and end of HP filtered, resampled data
  firstIndex = 1 + params.bufferSecs*params.resampleRate;
  lastIndex  = length(highpassed.data)- ...
    params.bufferSecs*params.resampleRate;

  r = constructTimeSeries(highpassed.data(firstIndex:lastIndex), ...
			  highpassed.tlow + params.bufferSecs, ...
			  highpassed.deltaT, ...
			  highpassed.fbase, highpassed.phase);

  % estimate power spectra for optimal filter
  [temp,freqs] = psd(r.data, params.PSD.FFTLength, ...
    1/r.deltaT, params.PSD.Window, params.PSD.OverlapLength, ...
    params.PSD.detrendFlag);

  % normalize appropriately: if all the bins in the PSD are independent, we
  % are dealing with complex params.heterodyned data
  spec = constructFreqSeries(2*r.deltaT*temp, freqs(1), ...
			     freqs(2)-freqs(1), 0);

  % coarse-grain noise power spectra to desired freqs
  PSD = coarseGrain(spec, params.flow, params.deltaF, params.numFreqs);

  % already calibrated, create FreqSeries structs
  calPSD = constructFreqSeries(PSD.data, ...
                            PSD.flow, PSD.deltaF, PSD.symmetry);

  params.midSegment = (params.numSegmentsPerInterval+1)/2;
  rbartilde = windowAndFFT(r, params.fft.dataWindow, ...
			       params.fft.fftLength);

return
