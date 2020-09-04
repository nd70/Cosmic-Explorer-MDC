function stochastic(paramsFile, startTime, jobDuration)
%
%  stochastic --- main routine for running the stochastic search
%
%  stochastic(paramsFile, startTime, jobDuration) executes the stochastic
%  search specified by a set of parameters on the chosen job.
%
%  As an alternative, for backwards compatability, if a jobs file is
%  available the pipeline can be executed with the parameters
%
%  stochastic(paramsFile, startTime, jobDuration)
%
%  - paramsFile is a text file containing one line for each parameter
%  to be set. Each line must be of the form
%
%    <parameters_name> <value>
%
%  The matlab comment symbol % may be used to include comments. Any
%  parameters not set will retain their default value. The file may
%  be empty or a null value for the file may be specified using the
%  value [], in which case all parameters will retain their defaults.
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
%  - startTime is the start time of the interval to be analyzed, as a
%    GPS time
%
%  - jobDuration is the duration of the interval in seconds. If the
%    duration is zero then the pipeline will exit immediately after
%    checking parameters without analyzing any data
%
%  - jobsFile is a text file containing one line for each job with the format
%
%    <jobNumber> <GPS_start_time> <GPS_stop_time> <duration>
%
%  Lines starting with the matlab comment symbol % are ignored.
%
%  - jobNumber must be an integer from 0 up to the number of non-commented
%  lines in jobsFile. A jobNumber of 0 is a "dummy" job which can be used to
%  quickly run the compiled executable to and check that the parameters are valid.
%
%  Routine written by Joseph D. Romano, John T. Whelan, Martin McHugh,
%  and Vuk Mandic.
%  Contact Joseph.Romano@astro.cf.ac.uk, john.whelan@ligo.org,
%  mmchugh@loyno.edu, and/or vmandic@ligo.caltech.edu
%  ethrane@physics.umn.edu for intermediate frame issues
%
%  $Id: stochastic.m,v 1.120 2009/02/11 00:29:20 ethrane Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check that the number of input parameters is correct
if (nargin ~= 3)
  fprintf('Usage: stochastic(paramsFile, startTime, jobDuration) or stochastic(paramsFile, jobsFile, jobNumber)\n');
  return;
end;

% Set the stochastic parameters. Parameters are set by default
% and then may be overridden by user-defined values in the params file.
params = setStochasticParams(paramsFile, startTime, jobDuration);

% Check the parameters and exit if errors are found
checkParamsStochastic(params);

% If params.jobDuration is zero, this is a dry run to verify parameters so it is
% not an error
if (params.jobDuration == 0)
  warning('Dummy job, all parameters checked - exiting.');
  return;
end;

% Set up various derived parameters (like full paths to the frame cache)
params = loadAuxiliaryInput(params);

% Check the frame cache, creating it if it is not found and useDatafindServer is true.
checkFrameCache(params);

% Set up spherical harmonics if necessary
if params.doSphericalHarmonics
  [w1w2bar, w1w2squaredbar, w1w2ovlsquaredbar] = windowFactors(params.fft1.dataWindow, params.fft2.dataWindow);
  vSpH = SpH(params.gammaLM_coeffsPath, [params.ifo1(1), params.ifo2(1) ], ...
           params.SpHLmax, params.numFreqs, params.flow,params.deltaF, ...
           params.Hf.data, params.mask.data, w1w2bar, w1w2squaredbar, ...
	   params.segmentDuration, params.SpHFreqIntFlag, ...
	   params.outputFilePrefix, ...
	   params.maxSegmentsPerMatfile, params.jobNumber);
  % dumpSpH(vSpH);
end

% Initialize live combining of results of isotropic analysis
initCombine(params);

% If we aren't using matlab output, write parameters to a text file
if (~params.writeOutputToMatFile)
  printParams(params);
end;

% Open files for saving results
params = openOutputFiles(params);

% simulate stochastic signals for the duration of the whole job
if params.doMonteCarlo
  initMonteCarlo(params);
end

isFirstPass = true;
% analyse the data
for I = 1:params.numIntervalsTotal

  badSegmentData = false;
  badResponse = false;

  if params.intermediate
    params.intervalStartTime = params.centeredStartTime + ...
      (I+1)*params.intervalTimeStride;
  else
    params.intervalStartTime = params.centeredStartTime + ...
      (I-1)*params.intervalTimeStride;
  end

  if ~params.intermediate %cet-------------------------------------------------
    % check if first pass through the loop
    if isFirstPass
    
      for J=1:params.numSegmentsPerInterval
    
        % read in time-series data from frames
        dataStartTime1 = params.intervalStartTime + (J-1)*params.segmentDuration - params.bufferSecs1;
        dataStartTime2 = params.intervalStartTime + (J-1)*params.segmentDuration - params.bufferSecs2;

        dataDuration1 = params.segmentDuration + 2*params.bufferSecs1;
        dataDuration2 = params.segmentDuration + 2*params.bufferSecs2;

        if dataStartTime1+dataDuration1 > params.lastLoadedDataEnd1
          params.lastLoadedDataEnd1 = min(dataStartTime1+params.minDataLoadLength,params.startTime+params.jobDuration);
          lastLoadedDataStart1 = dataStartTime1;
          tmpDuration = params.lastLoadedDataEnd1 - lastLoadedDataStart1;
          [longadcdata1, data1OK] = readTimeSeriesData(params.channelName1, ...
                                                       dataStartTime1, tmpDuration, ...
                                                       params.frameDuration1, ...
                                                       params.gpsTimesFile1, params.frameCacheFile1);
          if params.doInjFromFile1
            [longinjdata1, injData1OK] = readTimeSeriesData(params.injChannelName1, ...
                                                            dataStartTime1, tmpDuration, ...
                                                            params.injFrameDuration1, ...
                                                            params.injGPSTimesFile1, params.injFrameCacheFile1);
            if ( longadcdata1.tlow ~= longinjdata1.tlow )
              error('injected and ADC data have different start times');
            end
            if ( longadcdata1.deltaT ~= longinjdata1.deltaT )
              error('injected and ADC data have different sampling rates');
            end
            longadcdata1.data = longadcdata1.data ...
                + params.injScale1 * longinjdata1.data;
          else
            injData1OK = true;
          end
        end

        if dataStartTime2+dataDuration2 > params.lastLoadedDataEnd2
          params.lastLoadedDataEnd2 = min(dataStartTime2+params.minDataLoadLength,params.startTime+params.jobDuration);
          lastLoadedDataStart2 = dataStartTime2;
          tmpDuration = params.lastLoadedDataEnd2 - lastLoadedDataStart2;
          [longadcdata2, data2OK] = readTimeSeriesData(params.channelName2, ...
                                                       dataStartTime2,tmpDuration, ...
                                                       params.frameDuration2, ...
                                                       params.gpsTimesFile2, params.frameCacheFile2);
          if params.doInjFromFile2
            [longinjdata2, injData2OK] = readTimeSeriesData(params.injChannelName2, ...
                                                            dataStartTime2, tmpDuration, ...
                                                            params.injFrameDuration2, ...
                                                            params.injGPSTimesFile2, params.injFrameCacheFile2);
            if ( longadcdata2.tlow ~= longinjdata2.tlow )
              error('injected and ADC data have different start times');
            end
            if ( longadcdata2.deltaT ~= longinjdata2.deltaT )
              error('injected and ADC data have different sampling rates');
            end
            longadcdata2.data = longadcdata2.data + ...
                + params.injScale2 * longinjdata2.data;
          else
            injData2OK = true;
          end
        end

        % if either data stream is bad, set flag and exit loop
        % this needs to be done before we start referring to the
        % variables longadcdata1 or longadcdata2
        if ( (data1OK==false) | (data2OK==false) | ...
             (injData1OK==false) | (injData2OK==false) )
           badSegmentData = true;
           break
        end

        startindex = (dataStartTime1 - lastLoadedDataStart1) ...
                      /longadcdata1.deltaT + 1;
        endindex = (dataStartTime1 + dataDuration1 - lastLoadedDataStart1) ...
                      /longadcdata1.deltaT;
        adcdata1 = longadcdata1;
        adcdata1.data = longadcdata1.data(startindex:endindex);
        adcdata1.tlow = dataStartTime1;

        startindex = (dataStartTime2 - lastLoadedDataStart2) ...
                      /longadcdata2.deltaT + 1;
        endindex = (dataStartTime2 + dataDuration2 - lastLoadedDataStart2) ...
                      /longadcdata2.deltaT;
        adcdata2 = longadcdata2;
        adcdata2.data = longadcdata2.data(startindex:endindex);
        adcdata2.tlow = dataStartTime2;

        % KLUDGE: can override base frequency in parameter file
        if (~isnan(params.fbase1) )
	  adcdata1.fbase = params.fbase1;
        end;
        if (~isnan(params.fbase2) )
	  adcdata2.fbase = params.fbase2;
        end;
        % End KLUDGE

        if ( isnan(adcdata1.fbase) & isnan(adcdata2.fbase) )
	  if params.heterodyned
	    error('Trying to do params.heterodyned analysis on non-params.heterodyned data');
          end
        else
	  if (~params.heterodyned)
	    error('Trying to do non-params.heterodyned analysis on params.heterodyned data');
  	  end
        end

        % apply timing glitch correction
        if params.doTimingTransientSubtraction1
          adcdata1.data=periodicSubtract(adcdata1.data,params.TimingTransient1.TimingTransient.data);
        end
        if params.doTimingTransientSubtraction2
          adcdata2.data=periodicSubtract(adcdata2.data,params.TimingTransient2.TimingTransient.data);
        end
      
        % downsample the data 
        n1(J) = downsampleTimeSeries(adcdata1, ...
                                     params.resampleRate1, ...
                                     params.nResample1, ...
                                     params.betaParam1);
        
        n2(J) = downsampleTimeSeries(adcdata2, ...
                                     params.resampleRate2, ...
                                     params.nResample2, ...
                                     params.betaParam2);
        
        % free-up some memory
        clear adcdata1; 
        clear adcdata2;

        % calculate response functions from calibration data
        calibsec1 = dataStartTime1 + params.bufferSecs1;
        [response1(J), responseOK1]=calculateResponseFunction(params.cal1, calibsec1, ...
                                                          params.flow, params.deltaF, ...
                                                          params.numFreqs, ...
                                                          params.ASQchannel1, ...
                                                          params.channel1Calibrated);

        calibsec2 = dataStartTime2 + params.bufferSecs2;
        [response2(J), responseOK2]=calculateResponseFunction(params.cal2, calibsec2, ...
                                                          params.flow, params.deltaF, ...
                                                          params.numFreqs, ...
                                                          params.ASQchannel2, ...
                                                          params.channel2Calibrated);

        if (~responseOK1 | ~responseOK2)
          badResponse = true;
          break
        end

      end % loop over segments J

      % if bad data or bad response function for any segment, continue with 
      % next interval
      if (badSegmentData | badResponse)
        continue
      else
        isFirstPass = false;
      end
    else % not first pass

      % shift data and response functions accordingly
      for J=1:params.numSegmentsPerInterval-1

        if params.doOverlap
          % shift data by half a segment; need to worry about buffer
          N1 = length(n1(J).data);
          N2 = length(n2(J).data);
          bufferOffset1 = params.bufferSecs1/n1(J).deltaT;
          bufferOffset2 = params.bufferSecs2/n2(J).deltaT;

          data  = [n1(J).data(N1/2+1-bufferOffset1:N1-bufferOffset1) ; ...
                   n1(J+1).data(1+bufferOffset1:N1/2+bufferOffset1)];
          tlow  = n1(J).tlow+params.intervalTimeStride;
          n1(J) = constructTimeSeries(data, tlow, n1(J).deltaT, ...
                                      n1(J).fbase, n1(J).phase);

          data  = [n2(J).data(N2/2+1-bufferOffset2:N2-bufferOffset2) ; ...
                   n2(J+1).data(1+bufferOffset2:N2/2+bufferOffset2)];
          tlow  = n2(J).tlow+params.intervalTimeStride;
          n2(J) = constructTimeSeries(data, tlow, n2(J).deltaT, ...
                                      n2(J).fbase, n2(J).phase);
          
          % get response function corresponding to shifted start times
          calibsec1 = n1(J).tlow + params.bufferSecs1;
          [response1(J), responseOK1]=calculateResponseFunction(params.cal1, calibsec1, ...
                                                            params.flow, params.deltaF, ...
                                                            params.numFreqs, ...
                                                            params.ASQchannel1, ...
                                                            params.channel1Calibrated);

          calibsec2 = n2(J).tlow + params.bufferSecs2;
          [response2(J), responseOK2]=calculateResponseFunction(params.cal2, calibsec2, ...
                                                            params.flow, params.deltaF, ...
                                                            params.numFreqs, ...
                                                            params.ASQchannel2, ...
                                                            params.channel2Calibrated);
          if (~responseOK1 | ~responseOK2)
            badResponse = true;
            break
          end

        else
          % simple shift by a full segment
          n1(J)=n1(J+1);
          n2(J)=n2(J+1);
          response1(J)=response1(J+1);
          response2(J)=response2(J+1);
        end
      end % loop over J
      
      % if bad response function for any segment, continue with next interval
      if badResponse
        continue
      end
  
      % read in time-series data for next segment
      dataStartTime1 = params.intervalStartTime ...
          + (params.numSegmentsPerInterval-1)*params.segmentDuration ...
          - params.bufferSecs1;
      dataStartTime2 = params.intervalStartTime ...
          + (params.numSegmentsPerInterval-1)*params.segmentDuration ...
          - params.bufferSecs2;

      dataDuration1 = params.segmentDuration + 2*params.bufferSecs1;
      dataDuration2 = params.segmentDuration + 2*params.bufferSecs2;

      if dataStartTime1+dataDuration1 > params.lastLoadedDataEnd1
        params.lastLoadedDataEnd1 = min(dataStartTime1+params.minDataLoadLength,params.startTime+params.jobDuration);
        lastLoadedDataStart1 = dataStartTime1;
        tmpDuration = params.lastLoadedDataEnd1 - lastLoadedDataStart1;
        [longadcdata1, data1OK] = readTimeSeriesData(params.channelName1, ...
                                                     dataStartTime1, tmpDuration, ...
                                                     params.frameDuration1, ...
                                                     params.gpsTimesFile1, params.frameCacheFile1);
        if params.doInjFromFile1
          [longinjdata1, injData1OK] = readTimeSeriesData(params.injChannelName1, ...
                                                          dataStartTime1, tmpDuration, ...
                                                          params.injFrameDuration1, ...
                                                          params.injGPSTimesFile1, params.injFrameCacheFile1);
          if ( longadcdata1.tlow ~= longinjdata1.tlow )
            error('injected and ADC data have different start times');
          end
          if ( longadcdata1.deltaT ~= longinjdata1.deltaT )
            error('injected and ADC data have different sampling rates');
          end
          longadcdata1.data = longadcdata1.data + ...
              + params.injScale1 * longinjdata1.data;
        else
          injData1OK = true;
        end
      end

      if dataStartTime2+dataDuration2 > params.lastLoadedDataEnd2
        params.lastLoadedDataEnd2 = min(dataStartTime2+params.minDataLoadLength,params.startTime+params.jobDuration);
        lastLoadedDataStart2 = dataStartTime2;
        tmpDuration = params.lastLoadedDataEnd2 - lastLoadedDataStart2;
        [longadcdata2, data2OK] = readTimeSeriesData(params.channelName2, ...
                                                     dataStartTime2, tmpDuration, ...
                                                     params.frameDuration2, ...
                                                     params.gpsTimesFile2, params.frameCacheFile2);
        if params.doInjFromFile2
          [longinjdata2, injData2OK] = readTimeSeriesData(params.injChannelName2, ...
                                                          dataStartTime2, tmpDuration, ...
                                                          params.injFrameDuration2, ...
                                                          params.injGPSTimesFile2, params.injFrameCacheFile2);
          if ( longadcdata2.tlow ~= longinjdata2.tlow )
            error('injected and ADC data have different start times');
          end
          if ( longadcdata2.deltaT ~= longinjdata2.deltaT )
            error('injected and ADC data have different sampling rates');
          end
          longadcdata2.data = longadcdata2.data + ...
              + params.injScale2 * longinjdata2.data;
        else
          injData2OK = true;
        end
      end

      % if either data stream is bad, set flag and exit loop
      % this needs to be done before we start referring to the
      % variables longadcdata1 or longadcdata2
      if ( (data1OK==false) | (data2OK==false) | ...
           (injData1OK==false) | (injData2OK==false) )
        badSegmentData = true;
        break
      end
      
      startindex = (dataStartTime1 - lastLoadedDataStart1) ...
          /longadcdata1.deltaT + 1;
      endindex = (dataStartTime1 + dataDuration1 - lastLoadedDataStart1) ...
          /longadcdata1.deltaT;
      adcdata1 = longadcdata1;
      adcdata1.data = longadcdata1.data(startindex:endindex);
      adcdata1.tlow = dataStartTime1;
      startindex = (dataStartTime2 - lastLoadedDataStart2) ...
          /longadcdata2.deltaT + 1;
      endindex = (dataStartTime2 + dataDuration2 - lastLoadedDataStart2) ...
          /longadcdata2.deltaT;
      adcdata2 = longadcdata2;
      adcdata2.data = longadcdata2.data(startindex:endindex);
      adcdata2.tlow = dataStartTime2;

      % KLUDGE: can override base frequency in parameter file
      if (~isnan(params.fbase1) )
        adcdata1.fbase = params.fbase1;
      end
      if (~isnan(params.fbase2) )
        adcdata2.fbase = params.fbase2;
      end
      % End KLUDGE
  
      if ( isnan(adcdata1.fbase) & isnan(adcdata2.fbase) )
        if params.heterodyned
          error('Trying to do params.heterodyned analysis on non-params.heterodyned data');
        end
      else
        if (~params.heterodyned)
          error('Trying to do non-params.heterodyned analysis on params.heterodyned data');
        end
      end
                                                                            
      % apply timing glitch correction
      if params.doTimingTransientSubtraction1
        adcdata1.data=periodicSubtract(adcdata1.data,params.TimingTransient1.TimingTransient.data);
      end
      if params.doTimingTransientSubtraction2
        adcdata2.data=periodicSubtract(adcdata2.data,params.TimingTransient2.TimingTransient.data);
      end
    
      % downsample the data 
      n1(params.numSegmentsPerInterval) = downsampleTimeSeries(adcdata1, ...
                                                        params.resampleRate1, ...
                                                        params.nResample1, ...
                                                        params.betaParam1);
      
      n2(params.numSegmentsPerInterval) = downsampleTimeSeries(adcdata2, ...
                                                        params.resampleRate2, ...
                                                        params.nResample2, ...
                                                        params.betaParam2);


      % free-up some memory
      clear adcdata1; 
      clear adcdata2;

      calibsec1 = dataStartTime1 + params.bufferSecs1;
      [response1(J), responseOK1]=calculateResponseFunction(params.cal1, calibsec1, ...
                                                        params.flow, params.deltaF, ...
                                                        params.numFreqs, ...
                                                        params.ASQchannel1, ...
                                                        params.channel1Calibrated);
      
      calibsec2 = dataStartTime2 + params.bufferSecs2;
      [response2(J), responseOK2]=calculateResponseFunction(params.cal2, calibsec2, ...
                                                        params.flow, params.deltaF, ...
                                                        params.numFreqs, ...
                                                        params.ASQchannel2, ...
                                                        params.channel2Calibrated);

      if (~responseOK1 | ~responseOK2)
        badResponse = true;
        break
      end

    end % of if isFirstPass ... else ... end

    % if bad data or bad response function for any segment, continue with 
    % next interval
    if (badSegmentData | badResponse)
      continue
    end

  end %cet --> closes  if ~params.intermediate

  % loop over number of trials: params.numTrials=1 if params.doMonteCarlo=false
  for K=1:params.numTrials
   if ~params.intermediate %cet------------------------------------------------
                           % initialize data array for average psds
   
 
      if params.doMedianPSD & ~params.ignoreMidSegment
         med_data1 = zeros(params.numFreqs,params.numSegmentsPerInterval);
         med_data2 = zeros(params.numFreqs,params.numSegmentsPerInterval);
      elseif params.doMedianPSD & params.ignoreMidSegment
         med_data1 = zeros(params.numFreqs,params.numSegmentsPerInterval-1);
         med_data2 = zeros(params.numFreqs,params.numSegmentsPerInterval-1);
      else
         avg_data1 = zeros(params.numFreqs,1);
         avg_data2 = zeros(params.numFreqs,1);
      end

      % loop over number of segments
      for J=1:params.numSegmentsPerInterval

        % time-shift the data
        if params.doShift1
          shiftoffset = round(params.ShiftTime1 / n1(J).deltaT);
          qtempdata1 = circshift(n1(J).data,shiftoffset);
        else
          qtempdata1 = n1(J).data;
        end
        if params.doShift2
          shiftoffset = round(params.ShiftTime2 / n2(J).deltaT);
          qtempdata2 = circshift(n2(J).data,shiftoffset);
        else
          qtempdata2 = n2(J).data;
        end
      
        o1 = n1(J);
        o2 = n2(J);
        o1.data = qtempdata1;
        o2.data = qtempdata2;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % simulate detector noise
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if params.doSimulatedDetectorNoise
          %fprintf('simulating detector noise\n');
          [sdn1, sdn2] = getDetectorNoiseData(o1.tlow,o1.deltaT*length(o1.data));
          if params.simulatedDetectorNoiseNoRealData
            o1 = constructTimeSeries(sdn1, ...
                                     o1.tlow,  o1.deltaT, ...
                                     o1.fbase, o1.phase);
            o2 = constructTimeSeries(sdn2, ...
                                     o2.tlow,  o2.deltaT, ...
                                     o2.fbase, o2.phase);
          else
            o1 = constructTimeSeries(o1.data + sdn1, ...
                                     o1.tlow,  o1.deltaT, ...
                                     o1.fbase, o1.phase);
            o2 = constructTimeSeries(o2.data + sdn2, ...
                                     o2.tlow,  o2.deltaT, ...
                                     o2.fbase, o2.phase);
          end;
        end;
      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inject simulated isotropic signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if params.doMonteCarlo
          [o1,o2]=getDataMonteCarlo(params,o1,o2,I,J,K);
        end %if params.doMonteCarlo


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inject PointSource signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if params.doSimulatedPointSource
          %fprintf(' Get Data from %d to %d\n',o1.tlow,o1.tlow+o1.deltaT*length(o1.data));
          [sps1, sps2] = getPointSourceData(o1.tlow,o1.deltaT*length(o1.data));
          if params.simulatedPointSourcesNoRealData
            o1 = constructTimeSeries(sps1, ...
                                     o1.tlow,  o1.deltaT, ...
                                     o1.fbase, o1.phase);
            o2 = constructTimeSeries(sps2, ...
                                     o2.tlow,  o2.deltaT, ...
                                     o2.fbase, o2.phase);
          else
            o1 = constructTimeSeries(o1.data + sps1, ...
                                     o1.tlow,  o1.deltaT, ...
                                     o1.fbase, o1.phase);
            o2 = constructTimeSeries(o2.data + sps2, ...
                                     o2.tlow,  o2.deltaT, ...
                                     o2.fbase, o2.phase);
	  end;
        end;
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % inject SkyMap signal (spherical harmonics in time domain)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if params.doSimulatedSkyMap && params.simulatedSkyMapInjectTimeDomain
          %fprintf(' Get Data from %d to %d\n',o1.tlow,o1.tlow+o1.deltaT*length(o1.data));
          [ssm1, ssm2] = getSkyMapData(o1.tlow,o1.deltaT*length(o1.data));
	  if params.simulatedSkyMapNoRealData
            o1 = constructTimeSeries(ssm1, ...
                                     o1.tlow,  o1.deltaT, ...
                                     o1.fbase, o1.phase);
            o2 = constructTimeSeries(ssm2, ...
                                     o2.tlow,  o2.deltaT, ...
                                     o2.fbase, o2.phase);
          else
            o1 = constructTimeSeries(o1.data + ssm1, ...
                                     o1.tlow,  o1.deltaT, ...
                                     o1.fbase, o1.phase);
            o2 = constructTimeSeries(o2.data + ssm2, ...
                                     o2.tlow,  o2.deltaT, ...
                                     o2.fbase, o2.phase);
          end;
        end;
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % high-pass filter the data in-place. If no filtering is performed the
        % data is unchanged
        if params.doHighPass1
          if (params.useCascadeFilter1)
            o1.data = sos_filter(o1.data, params.highPassOrder1, ...
                                 params.highPassFreq1, params.resampleRate1);
          else
            o1.data = filtfilt(params.filt1.b, params.filt1.a, o1.data);
          end
        end

        if params.doHighPass2
          if (params.useCascadeFilter2)
            o2.data = sos_filter(o2.data, params.highPassOrder2, ...
                                 params.highPassFreq2, params.resampleRate2);
          else
            o2.data = filtfilt(params.filt2.b, params.filt2.a, o2.data);
          end
        end

        % chop-off bad data at start and end of HP filtered, resampled data
        firstIndex1 = 1 + params.bufferSecs1*params.resampleRate1;
        firstIndex2 = 1 + params.bufferSecs2*params.resampleRate2;

        lastIndex1  = length(o1.data) - params.bufferSecs1*params.resampleRate1;
        lastIndex2  = length(o2.data) - params.bufferSecs2*params.resampleRate2;

        r1(J) = constructTimeSeries(o1.data(firstIndex1:lastIndex1), ...
                                    o1.tlow + params.bufferSecs1, ...
                                    o1.deltaT, ...
                                    o1.fbase, o1.phase);
        r2(J) = constructTimeSeries(o2.data(firstIndex2:lastIndex2), ...
                                    o2.tlow + params.bufferSecs2, ...
                                    o2.deltaT, ...
                                    o2.fbase, o2.phase);

        % estimate power spectra for optimal filter
        calPSD1 = estimatePowerSpectrum(r1(J), params.psd1.FFTLength, ...
                                        params.psd1.Window, ...
                                        params.psd1.OverlapLength, ...
                                        response1(J));
        
        calPSD2 = estimatePowerSpectrum(r2(J), params.psd2.FFTLength, ...
                                        params.psd2.Window, ...
                                        params.psd2.OverlapLength, ...
                                        response2(J));


        if params.doMedianPSD
           % stores psds in order to calculate median power spectra, ignoring middle segment if desired
           params.midSegment = (params.numSegmentsPerInterval+1)/2;
          if ( (params.ignoreMidSegment) & (J==params.midSegment) )
            % do nothing
            %fprintf('Ignoring middle segment\n');
          else
           med_data1(:,J)=calPSD1.data;
           med_data2(:,J)=calPSD2.data;
          end
        else 
          % calculate avg power spectra, ignoring middle segment if desired
          params.midSegment = (params.numSegmentsPerInterval+1)/2;
          if ( (params.ignoreMidSegment) & (J==params.midSegment) )
            % do nothing
            %fprintf('Ignoring middle segment\n');
          else
            avg_data1 = avg_data1 + calPSD1.data;
            avg_data2 = avg_data2 + calPSD2.data;
          end
        end % end of if params.doMedianPSD

          if ( (params.writeNaiveSigmasToFiles) & (J==params.midSegment) )
	  % This calculates the "naive" theorerical variance, i.e.,
   	  % that calculated from the current segment without averaging 
	  % over the whole interval.
	  % This is useful for the stationarity veto which excludes
          % segments for which the naive sigma differs too much from
	  % the one calculated with the sliding PSD average.
	  [naiQ, result.naiVar, naiSensInt] ...
	      = calOptimalFilter(params.segmentDuration, params.gamma, ...
                                 params.fRef, params.alphaExp, ... 
                                 calPSD1, calPSD2, ...
                                 params.fft1.dataWindow, params.fft2.dataWindow, params.mask);
          if params.doAllSkyComparison
	      % Also do the same thing for the whole sky
	      [naiQAllSky, result.naiVarAllSky, naiSensIntAllSky] ...
	          = calOptimalFilter(params.segmentDuration, params.gammaAllSky, ...
                                     params.fRef, params.alphaExp, ... 
                                     calPSD1, calPSD2, ...
                                     params.fft1.dataWindow, params.fft2.dataWindow, params.mask);
  	    end
          end
 
        if ( (params.writeCoherenceToFiles) & (J==params.midSegment) )
          numOverlapInt= 2*params.numIntervalsTotal - 1;
          halfsize1=floor(length(r1(J).data)/2);
          halfsize2=floor(length(r2(J).data)/2);

          if I==1
            [c12,p1,p2, params.coh.f_coh] = cohereTimeSeries(r1(J).data, r2(J).data, ...
                                                             r1(J).fbase, r2(J).fbase, ...
                                                             r1(J).deltaT, r2(J).deltaT, ...
                                                             params.flow, params.fhigh);
            r1J_old=r1(J).data(halfsize1+1:length(r1(J).data));
            r2J_old=r2(J).data(halfsize2+1:length(r2(J).data));

            params.coh.cavg=c12/numOverlapInt;
            params.coh.p1avg=p1/numOverlapInt;
            params.coh.p2avg=p2/numOverlapInt;
            params.coh.coherenceStartTime = r1(J).tlow;
          else
            r1J_new = [r1J_old;r1(J).data(1:halfsize1)];
            r2J_new = [r2J_old;r2(J).data(1:halfsize2)];
            [c12,p1,p2, params.coh.f_coh] = cohereTimeSeries(r1J_new, r2J_new, ...
                                                             r1(J).fbase, r2(J).fbase, ...
                                                             r1(J).deltaT, r2(J).deltaT, ...
                                                             params.flow, params.fhigh);
            r1J_old=r1(J).data(halfsize1+1:length(r1(J).data));
            r2J_old=r2(J).data(halfsize2+1:length(r2(J).data));

            params.coh.cavg=params.coh.cavg + c12/numOverlapInt;
            params.coh.p1avg=params.coh.p1avg + p1/numOverlapInt;
            params.coh.p2avg=params.coh.p2avg + p2/numOverlapInt;
            [c12,p1,p2, params.coh.f_coh] = cohereTimeSeries(r1(J).data, r2(J).data, ...
                                                             r1(J).fbase, r2(J).fbase, ...
                                                             r1(J).deltaT, r2(J).deltaT, ...
                                                             params.flow, params.fhigh);
            params.coh.cavg=params.coh.cavg + c12/numOverlapInt;
            params.coh.p1avg=params.coh.p1avg + p1/numOverlapInt;
            params.coh.p2avg=params.coh.p2avg + p2/numOverlapInt;
          end
          
          params.coh.coh_avg=(abs(params.coh.cavg).^2)./(params.coh.p1avg.*params.coh.p2avg);     
        end
      
      end % loop over segments J


      if params.doMedianPSD
        % construct psd by taking median across times in each frequency bin
        if params.ignoreMidSegment
           % crop out middle segment
           idx_mid=(params.numSegmentsPerInterval-1)/2;
           med_data1=[med_data1(:,1:idx_mid-1),med_data1(:,idx_mid+1:end)];
           med_data2=[med_data2(:,1:idx_mid-1),med_data2(:,idx_mid+1:end)];
        else
                
        end
        % bias factor for using median -- see gr-qc/0509116 Appendix B

        % finite n -- only works for odd n
        %nsegs=length(med_data1(1,:));
        %ll_bias=1:nsegs;
        %alpha_bias=sum((-1).^(ll_bias+1)./ll_bias);

        % theoretical bias factor for infinite n
        alpha_bias=1; % FUDGE

        avg_data1 = median(med_data1,2) / alpha_bias;
        avg_data2 = median(med_data2,2) / alpha_bias;
      else
        % construct average power spectra
        if params.ignoreMidSegment
          avg_data1 = avg_data1/(params.numSegmentsPerInterval-1);
          avg_data2 = avg_data2/(params.numSegmentsPerInterval-1);
        else
          avg_data1 = avg_data1/params.numSegmentsPerInterval;
          avg_data2 = avg_data2/params.numSegmentsPerInterval;
        end
      end


 
      result.calPSD1_avg = constructFreqSeries(avg_data1, params.flow, params.deltaF, 0);
      result.calPSD2_avg = constructFreqSeries(avg_data2, params.flow, params.deltaF, 0);

    end %cet ---> if ~params.intermediate

    %CET: INTERMEDIATE FRAME STUFF HERE-------------------------------
    if params.intermediate
      %Define variables based on paramfile values.-----------------------------
      frameCacheFile =  params.intFrameFile;
      frameFiles = textread(frameCacheFile, '%s', -1, 'commentstyle', 'matlab');
      inputfile = char(frameFiles{I});
      ifo1=params.ifo1;
      ifo2=params.ifo2;
      gpsStart = params.intervalStartTime;
      segmentDuration = params.segmentDuration;
      params.midSegment = (params.numSegmentsPerInterval+1)/2;
      r1(params.midSegment).tlow = gpsStart;

      %EHT: read IM data from GPS=0
      gpsStart0 = 0;

      %Get intermediate data.--------------------------------------------------
      P1 = frgetvect(inputfile,[ifo1 ':AdjacentPSD'], ...
                     gpsStart0, segmentDuration);
      P2 = frgetvect(inputfile,[ifo2 ':AdjacentPSD'], ...
                     gpsStart0, segmentDuration);
      naiP1 = frgetvect(inputfile,[ifo1 ':LocalPSD'], ...
                        gpsStart0, segmentDuration);
      naiP2 = frgetvect(inputfile,[ifo2 ':LocalPSD'], ...
			gpsStart0, segmentDuration);
      CC = frgetvect(inputfile,[ifo1 ifo2 ':CSD'], ...
		     gpsStart0, segmentDuration);
      pp.flow = frgetvect(inputfile,[ifo1 ifo2 ':flow'], ...
                          gpsStart0, segmentDuration);
      pp.fhigh = frgetvect(inputfile,[ifo1 ifo2 ':fhigh'], ...
                           gpsStart0, segmentDuration);
      pp.deltaF = frgetvect(inputfile,[ifo1 ifo2 ':deltaF'], ...
                            gpsStart0, segmentDuration);
      pp.w1w2bar = frgetvect(inputfile,[ifo1 ifo2 ':w1w2bar'], ...
			     gpsStart0, segmentDuration);
      pp.w1w2squaredbar = frgetvect(inputfile, [ifo1 ifo2 ':w1w2squaredbar'], ...
				    gpsStart0, segmentDuration);
      pp.w1w2ovlsquaredbar = frgetvect(inputfile, [ifo1 ifo2 ':w1w2ovlsquaredbar'], ...
				      gpsStart0, segmentDuration);

      % May 5 - ethrane
      pp_freqs = pp.flow:pp.deltaF:pp.fhigh;
      pp_freq_cut = pp_freqs >= params.flow & pp_freqs <= params.fhigh;
      
      %      freqsToRemove = params.freqsToRemove;
      %      nBinsToRemove = params.nBinsToRemove;
      %      doFreqMask = params.doFreqMask;
      %      dataFM = constructFreqMask(pp.flow, pp.fhigh, pp.deltaF, ...
      %                         freqsToRemove, nBinsToRemove, doFreqMask);
      %      mask = constructFreqSeries(dataFM, pp.flow, pp.deltaF);
      %      numFreqs = floor((pp.fhigh-pp.flow)/pp.deltaF)+1;
      %      f = pp.flow + pp.deltaF*transpose([0:numFreqs-1]);
      %      dataORF = ...
      %	overlapreductionfunction(f, params.detector1, params.detector2);
      %      if params.heterodyned
      %        gamma = constructFreqSeries(dataORF, pp.flow, pp.deltaF, 0);
      %	else
      %        gamma = constructFreqSeries(dataORF, pp.flow, pp.deltaF, 1);
      %      end
      %      fRef = params.fRef;
      %      alphaExp = params.alphaExp;

      result.calPSD1_avg = constructFreqSeries(P1,pp.flow,pp.deltaF,1);
      result.calPSD2_avg = constructFreqSeries(P2,pp.flow,pp.deltaF,1);

      %cethrane: adding naiPSD_avg calcs for badgps work on August 3, 2009
      naiPSD1_avg = constructFreqSeries(naiP1,pp.flow,pp.deltaF,1);
      naiPSD2_avg = constructFreqSeries(naiP2,pp.flow,pp.deltaF,1);

      % May 5: ethrane - shorten the array for the case when
      % fhigh in processing does not match fhigh in pre-processing
      result.calPSD1_avg.data = result.calPSD1_avg.data(pp_freq_cut);
      result.calPSD2_avg.data = result.calPSD2_avg.data(pp_freq_cut);
      naiPSD1_avg.data = naiPSD1_avg.data(pp_freq_cut);
      naiPSD2_avg.data = naiPSD2_avg.data(pp_freq_cut);

      % ...and adjust flow to param value.
      result.calPSD1_avg.flow=params.flow;
      result.calPSD2_avg.flow=params.flow;
      naiPSD1_avg.flow=params.flow;
      naiPSD2_avg.flow=params.flow;

      %      [result.Q, result.ccVar, result.sensInt] = calOptimalFilter(segmentDuration, gamma, ...
      %					     fRef, alphaExp, ...
      %					     result.calPSD1_avg, result.calPSD2_avg, ...
      %					     -1, -1, mask, pp);
      %      maxCorrelationTimeShift=params.maxCorrelationTimeShift;
      %      CSD = constructFreqSeries(CC /2 * segmentDuration ...
      %				  * pp.w1w2bar,pp.flow,pp.deltaF,1);
      %      [result.ccStat,result.ccSpec] = processCSD(CSD, result.Q);
    end %cet --------------------------------------> if params.intermediate

    % calculate optimal filter, theoretical variance, and sensitivity 
    % integrand using avg psds
    [result.Q, result.ccVar, result.sensInt] = calOptimalFilter(params.segmentDuration, ...
                                                      params.gamma, ...
                                                      params.fRef, params.alphaExp, ... 
                                                      result.calPSD1_avg, result.calPSD2_avg, ...
                                                      params.fft1.dataWindow, ...
                                                      params.fft2.dataWindow, ...
                                                      params.mask);
    
    %Sept 15, 2009: take these lines b/c of error in S6 online studies
    %cethrane: adding naiVar calculation on August 3, 2009 for bad gps calcs.
    %     [naiQ, result.naiVar, naiSensInt] = calOptimalFilter(params.segmentDuration, ...
    %                                             params.gamma, ...
    %                                             params.fRef, params.alphaExp, ...
    %                                             naiPSD1_avg, naiPSD2_avg, ...
    %                                            params.fft1.dataWindow, ...
    %                                             params.fft2.dataWindow, ...
    %                                             params.mask);
    
    if params.doAllSkyComparison
      % Also do the same thing for the whole sky
      [QAllSky, result.ccVarAllSky, result.sensIntAllSky] = ...
          calOptimalFilter(params.segmentDuration, ...
                           params.gammaAllSky, ...
                           params.fRef, params.alphaExp, ... 
                           result.calPSD1_avg, result.calPSD2_avg, ...
                           params.fft1.dataWindow, ...
                           params.fft2.dataWindow, ...
                           params.mask);
    end
    if params.intermediate %cet
      CSD = constructFreqSeries(CC /2 * segmentDuration * pp.w1w2bar, ...
                                pp.flow,pp.deltaF, 1);

      % May 5: ethrane - shorten the array for the case when
      % fhigh in processing does not match fhigh in pre-processing
      CSD.data = CSD.data(pp_freq_cut);
      CSD.flow = params.flow;
      
      if params.doDirectional
        [result.ccStat,result.ccSpec] = processCSD(CSD, result.Q, params.maxCorrelationTimeShift);
      else
        [result.ccStat,result.ccSpec] = processCSD(CSD, result.Q);
      end
      if params.doAllSkyComparison
        [result.ccStatAllSky,result.ccSpecAllSky] = processCSD(CSD, result.QAllSky);
      end
    end
    
    %cet
    if ~params.intermediate
      % analyse the middle data segment
      % window, zero-pad and fft the data
      rbartilde1 = windowAndFFT(r1(params.midSegment), ...
                                params.fft1.dataWindow, ...
                                params.fft1.fftLength);
      rbartilde2 = windowAndFFT(r2(params.midSegment), ...
                                params.fft2.dataWindow, ...
                                params.fft2.fftLength);
      result.ccStat=NaN;
      dummydata=zeros(params.numFreqs,1).*NaN;
      result.ccSpec= constructFreqSeries(dummydata, params.flow, params.deltaF);
      if params.doDirectional
        [result.ccStat,result.ccSpec] = calCrossCorr(rbartilde1, rbartilde2, result.Q, ...
                                                     response1(params.midSegment), ...
                                                     response2(params.midSegment), ...
                                                     params.maxCorrelationTimeShift);
      else
        [result.ccStat,result.ccSpec] = calCrossCorr(rbartilde1, rbartilde2, result.Q, ...
                                                     response1(params.midSegment), ...
                                                     response2(params.midSegment));
      end
      if params.doAllSkyComparison
        [result.ccStatAllSky,result.ccSpecAllSky] = calCrossCorr(rbartilde1, ...
                                                          rbartilde2, result.QAllSky, ...
                                                          response1(params.midSegment), ...
                                                          response2(params.midSegment));
      end
    end % if ~params.intermediate
    
    % calculate the value and spectrum of the CC statistic
    if params.doSphericalHarmonics
      %cet: how to do callSpH with SID
      if ~params.intermediate
        vSpH = callSpH(vSpH,params,rbartilde1,rbartilde2,...
                       result.calPSD1_avg,result.calPSD2_avg,r1(params.midSegment).tlow);
      end
      if params.intermediate
        % Feb 24: ethrane - the next two lines shorten the array for the case when
        % fhigh in processing does not match fhigh in pre-processing
        %          CC = CC(1:length(params.mask.data));
        %          P1 = P1(1:length(params.mask.data));
        %          P2 = P2(1:length(params.mask.data));
        % May 5: ethrane - the next two lines shorten the array for the case when
        % fhigh in processing does not match fhigh in pre-processing
        CC = CC(pp_freq_cut);
        P1 = P1(pp_freq_cut);
        P2 = P2(pp_freq_cut);
        % apr 12 (cethrane: can not shift IM data, so do this)
        %          if params.IMScramblePhase
        %   	    CC = exp(2*pi*i*rand(size(CC))) .* CC;
        %          end
        vSpH = doSpH(vSpH,params,CC,P1,P2,gpsStart);
      end
      % assign dummy for now
      %   	  result.ccStat=NaN;
      %	  dummydata=zeros(params.numFreqs,1).*NaN;
      %	  result.ccSpec= constructFreqSeries(dummydata, params.flow, params.deltaF);
      %          [result.ccStat,result.ccSpec] = calCrossCorr(rbartilde1, rbartilde2, result.Q, ...
      %                                     response1(params.midSegment), ...
      %                                     response2(params.midSegment));
    elseif params.doDirectional
      %        [result.ccStat,result.ccSpec] = calCrossCorr(rbartilde1, rbartilde2, result.Q, ...
      %                                     response1(params.midSegment), ...
      %                                     response2(params.midSegment), ...
      %				       params.maxCorrelationTimeShift);
      % if in radiometer mode: still have to read out the time series ccStat
      map.time=r1(params.midSegment).tlow;
      GPSmid=map.time+params.segmentDuration/2;
      if params.doNarrowbandRadiometer
        map.data=ccSpecReadout(params.detector1, params.detector2, ...
                               GPSmid, params.SkyPattern, result.ccSpec, result.ccVar, ...
                               result.sensInt, params.UnphysicalTimeShift);
        
        map.data([1:size(params.SkyPattern,1)]+size(map.data,1),:) ...
            = ccStatReadout(params.detector1, params.detector2, GPSmid, ...
                            params.SkyPattern, result.ccStat, result.ccVar, ...
                            params.UnphysicalTimeShift);
        
        map.isNarrowbandRadiometer=true;
        map.flow                  =result.ccSpec.flow;
        map.deltaF                =result.ccSpec.deltaF;
        map.numFreqs              =length(result.ccSpec.data);
      else
        map.data=ccStatReadout(params.detector1, params.detector2, ...
                               GPSmid, params.SkyPattern, result.ccStat, result.ccVar, ...
                               params.UnphysicalTimeShift);
      end
      if params.doAllSkyComparison
        %          [result.ccStatAllSky,result.ccSpecAllSky] = ...
        %                        calCrossCorr(rbartilde1, rbartilde2, result.QAllSky, ...
        %                                     response1(params.midSegment), ...
        %                                     response2(params.midSegment));
        map.data(1+size(map.data,1),:)=[result.ccStatAllSky,sqrt(result.ccVarAllSky)];
        map.lastLineIsAllSky=true;
      end
      metamap.time=map.time;
      metamap.filename=params.ccStatSkySetCurrentFilename;
      metamap.segmentOffset=params.skySetSegmentOffset;
      params.metaSky{params.skyIndex}=metamap;
      params.Sky{params.skyIndex-params.skySetSegmentOffset}=map;
      params.skyIndex=params.skyIndex+1;
      if params.skyIndex-params.skySetSegmentOffset > params.maxSegmentsPerMatfile
        Sky=params.Sky;
        save(params.ccStatSkySetCurrentFilename,'Sky');
        clear Sky;
        params.Sky={};
        params.skySetSegmentOffset=params.skySetSegmentOffset+...
	    params.maxSegmentsPerMatfile;
        params.skySetNumber=params.skySetNumber+1;
        params.ccStatSkySetCurrentFilename = [params.ccStatSkySetPrefix ...
                            num2str(params.skySetNumber) ...
                            params.ccStatSkySetSuffix];
      end % ---> params.doSphericalHarmonics = true
    else
      %          [result.ccStat,result.ccSpec] = calCrossCorr(rbartilde1, rbartilde2, result.Q, ...
      %                                     response1(params.midSegment), ...
      %                                     response2(params.midSegment));
    end
    
    % DISPLAY, WRITE RESULTS TO FILES (if desired)
    writeToOutputFiles(params,I,K,r1(params.midSegment).tlow, result);

    %SHIFT: have to worry about different trials
    addCombine(params,I,r1(params.midSegment).tlow, result);
    
  end % loop over trials K

end % loop over intervals I

%SHIFT: have to worry about different trials
K = 1; % warning: not implemented for more than 1 trial
finishCombine(params, K);
 
if params.doSphericalHarmonics
  vSpH = saveSpHSet(vSpH);
end

% Write final data common to all intervals and close files
writeToOutputFilesFinal(params, K)
params = closeOutputFiles(params);

return; % end of stochastic()

function [response, responseOK]=calculateResponseFunction(cal, calibsec, ...
                                                    flow, deltaF, numFreqs, ...
                                                    channelName, channelCalibrated)
%
% Calculate the response function at a particular time and return it as a 
% frequency-series struct
%
% Inputs:
%   cal - calibration data structure
%   calibsec - the time (in GPS seconds) at which to calculate the response function
%   flow - initial frequency (Hz)
%   deltaF - frequency spacing (Hz)
%   numFreqs - desired number of frequencies in the response function
%   channelName - the name of the channel. This is needed to calculate the
%     appropriate response function
%   channelCalibrated - if true, indicates that the data being processed is
%     already calibrated, so the function will generate a trivial (all 1's)
%     response function
%
% Outputs:
%   response - the response function as a frequency-series struct
%   responseOK - if true, indicates that the calibration at calibsec was
%     considered "good". A value of false for this flag
%     indicates that the alpha value at the required was too small (<10^(-6))
%     for the calibration to be considered "good". A response function of 0 is
%     returned in this instance.
%
  if (channelCalibrated)
    % the data is already calibrated
    response = constructFreqSeries(ones(numFreqs, 1), flow, deltaF);
    responseOK = true;
    transfer = response;
  else
    [R, responseOK] = calculateResponse(cal.t, cal.f, cal.R0, cal.C0, ...
                                        cal.alpha, cal.gamma, calibsec, channelName);

    if (responseOK)
      % evaluate response function at desired frequencies
      response = convertResponse(cal.f, R, flow, deltaF, numFreqs, 0, 0);
    else
      response = constructFreqSeries(R, flow, deltaF);
    end;
  end;

return;

function calPSD=estimatePowerSpectrum(tseries, FFTLength, Window, ...
                                      OverlapLength, response)
%
% Estimate the power spectrum of a time-series object, calibrated with the
% supplied response function.
%
% Inputs:
%   tseries - a time-series structure created with chanstruct()
%   FFTLength - an integer specifying the FFT length to use
%   Window - a vector of length FFTLength specifying the window to apply
%   OverlapLength - the overlap to use when estimating the PSD
%   response - a frequency-series structure containing the response function
%
% Outputs:
%   calPSD - a frequency-series object containing the PSD of the data in
%   tseries. The PSD has been calibrated by being multiplied by |response|^2
%

  [psdTemp, freqs] = pwelch(tseries.data, Window, OverlapLength, FFTLength, 1/tseries.deltaT);

  deltaF = freqs(2) - freqs(1);

  % Normalize appropriately.  If all the bins in the PSD are independent, we are 
  % dealing with complex params.heterodyned data
  if (length(psdTemp) == FFTLength)
    % Account for heterodyning of data.
    % Pwelch() produces a "correct" two-sided spectrum (unlike psd()) which is
    % half the one-sided spectrum except at DC and Nyquist. Since we previously
    % used psd(), the code dealing with heterodyning expects the spectrum to be
    % twice what pwelch() produces.
    freqs_shifted = fftshift(freqs);
    spec = constructFreqSeries(2*fftshift(psdTemp), ...
                               tseries.fbase + freqs_shifted(1) - 1/tseries.deltaT, ...
                               deltaF, 0);
  else
    % One-sided spectrum
    spec = constructFreqSeries(psdTemp, freqs(1), ...
                               deltaF, 0);
  end

  % coarse-grain noise power spectra to desired freqs - these must match the
  % properties of the response function
  calPSD = coarseGrain(spec, response.flow, response.deltaF, length(response.data));

  % calibrate the power spectrum by applying the response function to the spectrum
  % values - other fields do no need to be changed
  calPSD.data = calPSD.data.*(abs(response.data).^2);

return;

function tseries=downsampleTimeSeries(tseries, resampleRate, nResample, betaParam)
%
% Downsample the data in a time-series structure
%
% This function takes a time-series structure and downsamples it to the new
% sampling rate provided in resampleRate. The parameters passed to Matlab's
% resample function are provided in nResample and betaParam.
%
% The result is returned in a time-series structure. All the fields of the new time-series
% are cloned from the input, other than the data itself which has been downsampled.
% This means, for example, that the start-time of the resulting time-series is the same as
% that of the input time-series - the delay of the filter used in the downsampling is
% not used to modify the start-time.
%
% If resampleRate is the same as the sampling rate of tseries, the time-series is
% returned unchanged
%
% If resampleRate is greater than the sampling rate of tseries, the function gives
% an error.
%
  sampleRate = 1/tseries.deltaT;
  p = 1;  
  q = floor(sampleRate/resampleRate);

  % Resample the time series in-place - if the new sampling rate is
  % the same as the old there is no need to do anything
  if (resampleRate < sampleRate)
    tseries.data = resample(tseries.data, p, q, nResample, betaParam);
    tseries.deltaT = 1/resampleRate;
  elseif (resampleRate > sampleRate)
    error(sprintf('New sample rate (%f) is higher than old sampling rate (%).', ...
                  resampleRate, sampleRate));
  end;

return;
