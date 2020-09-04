function [h1, h2, badResponse] = getDetectorNoiseData(GPSstart,dur)
%  function [h1, h2] = getDetectorNoiseData(GPSstart,dur)
%
%  getDetectorNoiseData --- returns dur seconds of simulated detector noise starting at GPSstart.
%                           The function relies on the GLOBAL variable DETECTOR_NOISE_MEMORY,
%                           which is a struct and has to be initialized properly before the 1st call.
%
%  arguments: GPSstart   - GPS start time of the data
%             dur        - duration of the data in seconds
%
%  DETECTOR_NOISE_MEMORY:
%             fsample1   - sample frequency of the time series 1
%             fsample2   - sample frequency of the time series 2
%             nResample1 - order of matlab resample routine for time series 1
%             nResample2 - order of matlab resample routine for time series 2
%             betaParam1 - beta parameter for matlab resample routine for time series 1
%             betaParam2 - beta parameter for matlab resample routine for time series 2
%             P1,2       - power spectra (one-sided, strain^2/Hz) in detectors 1,2
%                          preferred input: Nx2 array (filename and freq. series
%                          work too ... in principle)
%             intLog     - boolean whether to interpolate P1,2 logarithmicly
%             flow             \
%             deltaF            > frequency information used for the main analysis
%             numFreqs         |
%             t1,f1,           \
%             R01, C01,         \
%             alpha1,            \
%             gamma1,             \
%             ASQchannel1          \
%             alphaBetaFile1        \
%             calCavGainFile1        \
%             calResponseFile1        \ Calibration information
%             t2,f2,                  |
%             R02, C02,              |
%             alpha2,               |
%             gamma2,              |
%             ASQchannel2         |
%             alphaBetaFile2     |
%             calCavGainFile2   |
%             calResponseFile2 |
%             bufstart    GPS start time of buffer
%             bufforget   memory duration of the buffer - this determines how far back
%                         the data can be recalled, i.e. the identical random series is produced
%             bufdur      number of seconds of data in the buffer - ready to be read out, <=bufforget
%             bufspline   duration in sec of 1/4 period of a cos/sin used to spline data
%                         2 x bufspline is also the time interal used to get new data
%                         from simulateSkyMap
%             buf1        the actual buffer for IFO1
%             buf2        the actual buffer for IFO2
%
%  output:    h1,h2      - time series seen by detector 1,2, calibrated in counts
%
%  Routine written by Stefan Ballmer, Joe Romano
%  Contact sballmer@ligo.mit.edu
%
%  $Id: getDetectorNoiseData.m,v 1.4 2008-09-25 17:30:52 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global DETECTOR_NOISE_MEMORY;
badResponse = false;

% is GPSstart too far back?
if GPSstart < DETECTOR_NOISE_MEMORY.bufstart
  error('cannot go that far into the past - increase bufforget');
end

% Do I need to get new data?
if GPSstart + dur > DETECTOR_NOISE_MEMORY.bufstart + DETECTOR_NOISE_MEMORY.bufdur
  % does the buffer overlap with the new past?
  if GPSstart+dur- DETECTOR_NOISE_MEMORY.bufforget < ...
                   DETECTOR_NOISE_MEMORY.bufstart + DETECTOR_NOISE_MEMORY.bufdur + DETECTOR_NOISE_MEMORY.bufspline
    % does the buffer get too long?
    nspline=ceil((GPSstart+dur-DETECTOR_NOISE_MEMORY.bufstart-DETECTOR_NOISE_MEMORY.bufdur)/DETECTOR_NOISE_MEMORY.bufspline);
    lost=DETECTOR_NOISE_MEMORY.bufdur + nspline*DETECTOR_NOISE_MEMORY.bufspline - DETECTOR_NOISE_MEMORY.bufforget;
    if lost > 0
      % crop the buffer
      indi1=1 + lost*DETECTOR_NOISE_MEMORY.fsample1;
      indi2=1 + lost*DETECTOR_NOISE_MEMORY.fsample2;
      DETECTOR_NOISE_MEMORY.bufstart = DETECTOR_NOISE_MEMORY.bufstart+lost;
      DETECTOR_NOISE_MEMORY.bufdur   = DETECTOR_NOISE_MEMORY.bufdur-lost;
      DETECTOR_NOISE_MEMORY.buf1     = DETECTOR_NOISE_MEMORY.buf1(indi1:end);
      DETECTOR_NOISE_MEMORY.buf2     = DETECTOR_NOISE_MEMORY.buf2(indi2:end);
      %fprintf('Buffer cropped\n');
    end
  else
    % the buffer is useless - delete it
    DETECTOR_NOISE_MEMORY.bufstart = GPSstart;
    DETECTOR_NOISE_MEMORY.bufdur   = -DETECTOR_NOISE_MEMORY.bufspline; % indicates that there is also no rolled iff data in buf
    DETECTOR_NOISE_MEMORY.buf1     = [];
    DETECTOR_NOISE_MEMORY.buf2     = [];
    %fprintf('Buffer erased\n');
  end
  % define sine and cosine arrays
  N1=2*DETECTOR_NOISE_MEMORY.bufspline*DETECTOR_NOISE_MEMORY.fsample1;
  N2=2*DETECTOR_NOISE_MEMORY.bufspline*DETECTOR_NOISE_MEMORY.fsample2;
  s1 = sin(pi*[0:N1/2-1]'/N1);
  s2 = sin(pi*[0:N2/2-1]'/N2);
  c1 = cos(pi*[0:N1/2-1]'/N1);
  c2 = cos(pi*[0:N2/2-1]'/N2);
  % keep filling the buffer until enough data is in the can
  while GPSstart + dur > DETECTOR_NOISE_MEMORY.bufstart + DETECTOR_NOISE_MEMORY.bufdur , 
    % start time of new chunk
    chunkstart  =DETECTOR_NOISE_MEMORY.bufstart+DETECTOR_NOISE_MEMORY.bufdur;
    calibsec    =chunkstart+DETECTOR_NOISE_MEMORY.bufspline;
    siderealtime=GPStoGreenwichMeanSiderealTime(calibsec);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% calibrate                                                                   %%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % calculate transfer functions appropriate for this and half of the 
      % next segment (need both of these response functions for splicing)
      if ( strncmp(DETECTOR_NOISE_MEMORY.alphaBetaFile1,   'none', length(DETECTOR_NOISE_MEMORY.alphaBetaFile1))   | ...
           strncmp(DETECTOR_NOISE_MEMORY.calCavGainFile1,  'none', length(DETECTOR_NOISE_MEMORY.calCavGainFile1))  | ...
           strncmp(DETECTOR_NOISE_MEMORY.calResponseFile1, 'none', length(DETECTOR_NOISE_MEMORY.calResponseFile1)) )
        % the data is already calibrated
        transfer1 = constructFreqSeries(ones(DETECTOR_NOISE_MEMORY.numFreqs,1), DETECTOR_NOISE_MEMORY.flow, DETECTOR_NOISE_MEMORY.deltaF);
      else
        [R1, responseOK1] = ...
          calculateResponse(DETECTOR_NOISE_MEMORY.t1,    DETECTOR_NOISE_MEMORY.f1,...
			    DETECTOR_NOISE_MEMORY.R01,   DETECTOR_NOISE_MEMORY.C01,...
			    DETECTOR_NOISE_MEMORY.alpha1,DETECTOR_NOISE_MEMORY.gamma1,...
			    calibsec,                  DETECTOR_NOISE_MEMORY.ASQchannel1);
        % if response function is bad, set flag and exit loop
        if (responseOK1 == false)
          badResponse = true; h1=[]; h2=[];
          fprintf('bad response for detector 1 in Point Source simulation\n');
          return; % exit function
        end  
        % convert to transfer functions (units: counts/strain) 
        transfer1 = convertResponse(DETECTOR_NOISE_MEMORY.f1,    R1, DETECTOR_NOISE_MEMORY.flow,...
				    DETECTOR_NOISE_MEMORY.deltaF,    DETECTOR_NOISE_MEMORY.numFreqs, 1, 0);
      end

      if ( strncmp(DETECTOR_NOISE_MEMORY.alphaBetaFile2,   'none', length(DETECTOR_NOISE_MEMORY.alphaBetaFile2))   | ...
           strncmp(DETECTOR_NOISE_MEMORY.calCavGainFile2,  'none', length(DETECTOR_NOISE_MEMORY.calCavGainFile2))  | ...
           strncmp(DETECTOR_NOISE_MEMORY.calResponseFile2, 'none', length(DETECTOR_NOISE_MEMORY.calResponseFile2)) )

        % the data is already calibrated
        transfer2 = constructFreqSeries(ones(DETECTOR_NOISE_MEMORY.numFreqs,1), DETECTOR_NOISE_MEMORY.flow, DETECTOR_NOISE_MEMORY.deltaF);
      else
        [R2, responseOK2] = ...
          calculateResponse(DETECTOR_NOISE_MEMORY.t2,    DETECTOR_NOISE_MEMORY.f2,...
			    DETECTOR_NOISE_MEMORY.R02,   DETECTOR_NOISE_MEMORY.C02,...
			    DETECTOR_NOISE_MEMORY.alpha2,DETECTOR_NOISE_MEMORY.gamma2,...
			    calibsec,             DETECTOR_NOISE_MEMORY.ASQchannel2);
        % if response function is bad, set flag and exit loop
        if (responseOK2 == false)
          badResponse = true; h1=[]; h2=[];
          fprintf('bad response for detector 2 in Point Source simulation\n');
          return; % exit function
        end  
        % convert to transfer functions (units: counts/strain) 
        transfer2 = convertResponse(DETECTOR_NOISE_MEMORY.f2,    R2, DETECTOR_NOISE_MEMORY.flow,...
				    DETECTOR_NOISE_MEMORY.deltaF,    DETECTOR_NOISE_MEMORY.numFreqs, 1, 0);
       end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% end of calibration                                                          %%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get new data;
    fprintf('Simulating detector noise between %d and %d\n',chunkstart,chunkstart+DETECTOR_NOISE_MEMORY.bufspline*2);
    [d1,d2] =  simulateDetectorNoise(DETECTOR_NOISE_MEMORY.bufspline*2,...
                                     GPSstart,...
                                     DETECTOR_NOISE_MEMORY.fsample1,...
                                     DETECTOR_NOISE_MEMORY.fsample2,...
                                     DETECTOR_NOISE_MEMORY.nResample1,...
                                     DETECTOR_NOISE_MEMORY.nResample2,...
                                     DETECTOR_NOISE_MEMORY.betaParam1,...
                                     DETECTOR_NOISE_MEMORY.betaParam2,...
			             DETECTOR_NOISE_MEMORY.P1,...
			             DETECTOR_NOISE_MEMORY.P2,...
                                     DETECTOR_NOISE_MEMORY.intLog,...
                                     transfer1,transfer2);
    if or(any(isnan(d1)),any(isnan(d2)))
      error('simulated detector noise data contains NaN');
    end
    % refill bufs with new data;
    if DETECTOR_NOISE_MEMORY.bufdur < 0
      %special case: buffer has been deleted
      indi1=1;
      indi2=1;
      indf1=N1/2;
      indf2=N2/2;
      DETECTOR_NOISE_MEMORY.buf1(indi1:indf1,1)= c1.*d1(N1/2+1:N1);
      DETECTOR_NOISE_MEMORY.buf2(indi2:indf2,1)= c2.*d2(N2/2+1:N2);
    else
      %regular case: still data in the buffer
      indi1=1 + DETECTOR_NOISE_MEMORY.bufdur * DETECTOR_NOISE_MEMORY.fsample1;
      indi2=1 + DETECTOR_NOISE_MEMORY.bufdur * DETECTOR_NOISE_MEMORY.fsample2;
      indf1=indi1 + N1/2 - 1;
      indf2=indi2 + N2/2 - 1;
      DETECTOR_NOISE_MEMORY.buf1(indi1:indf1,1)=DETECTOR_NOISE_MEMORY.buf1(indi1:indf1) + s1.*d1(1:N1/2);
      DETECTOR_NOISE_MEMORY.buf2(indi2:indf2,1)=DETECTOR_NOISE_MEMORY.buf2(indi2:indf2) + s2.*d2(1:N2/2);
      indi1=indi1+N1/2;
      indi2=indi2+N2/2;
      indf1=indf1+N1/2;
      indf2=indf2+N2/2;
      DETECTOR_NOISE_MEMORY.buf1(indi1:indf1,1)= c1.*d1(N1/2+1:N1);
      DETECTOR_NOISE_MEMORY.buf2(indi2:indf2,1)= c2.*d2(N2/2+1:N2);
    end;
    % update pointer
    DETECTOR_NOISE_MEMORY.bufdur=DETECTOR_NOISE_MEMORY.bufdur + DETECTOR_NOISE_MEMORY.bufspline;
  end;
end;
% requested data is now in the buffer, return it
indi1=1 + (GPSstart-DETECTOR_NOISE_MEMORY.bufstart)*DETECTOR_NOISE_MEMORY.fsample1;
indi2=1 + (GPSstart-DETECTOR_NOISE_MEMORY.bufstart)*DETECTOR_NOISE_MEMORY.fsample2;
indf1=indi1 + dur*DETECTOR_NOISE_MEMORY.fsample1 - 1;
indf2=indi2 + dur*DETECTOR_NOISE_MEMORY.fsample2 - 1;
h1=DETECTOR_NOISE_MEMORY.buf1(indi1:indf1);
h2=DETECTOR_NOISE_MEMORY.buf2(indi2:indf2);

return;

