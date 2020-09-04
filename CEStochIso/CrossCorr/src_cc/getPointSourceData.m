function [h1, h2,badResponse] = getPointSourceData(GPSstart,dur)
%  function [h1, h2] = getPointSourceData(GPSstart,dur)
%
%  getPointSourceData --- returns dur seconds of simulated point source data starting at GPSstart
%                         The function relies on the GLOBAL variable POINT_SOURCE_MEMORY,
%                         which is a struct and has to be initialized properly before the 1st call.
%
%  arguments: GPSstart   - GPSt start time of the data
%             dur        - duration of the data in seconds
%
%  POINT_SOURCE_MEMORY:
%             fsample    - sample frequency of the time series
%             Hf         - total power spectrum (one-sided) in both polarizations, i.e.
%                          the actual power spectrum for each polarization is Hf/2
%                          preferred input: Nx2 array (filename and freq. series
%                          work too ... in principle)
%             intLog     - boolean whether to interpolate Hf logarithmicly
%             det1,det2  - detector structures containing position r and tensor d
%             ra         - right ascension in hours of source in the sky
%                          takes vector for multiple point sources, must have same length as decl
%             decl       - right ascension in degrees of source in the sky
%                          takes vector for multiple point sources, must have same length as ra
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
%                         from simulatePointSource
%             buf1        the actual buffer for IFO1
%             buf2        the actual buffer for IFO2
%             MakeIncoherent optional parameter to destroy coherence of point source
%                            0:  coherent point source
%                            1:  incoherent, but scaled with DC antenna acceptance
%                            2:  stationary noise, PowerSpec corresponds to the noise seen in both 
%
%  output:    h1,h2      - time series seen by detector 1,2, calibrated in counts
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global POINT_SOURCE_MEMORY;
badResponse = false;

% is GPSstart too far back?
if GPSstart < POINT_SOURCE_MEMORY.bufstart
  error('cannot go that far into the past - increase bufforget');
end

% Do I need to get new data?
if GPSstart + dur > POINT_SOURCE_MEMORY.bufstart + POINT_SOURCE_MEMORY.bufdur
  % does the buffer overlap with the new past?
  if GPSstart+dur- POINT_SOURCE_MEMORY.bufforget < ...
                   POINT_SOURCE_MEMORY.bufstart + POINT_SOURCE_MEMORY.bufdur + POINT_SOURCE_MEMORY.bufspline
    % does the buffer get too long?
    nspline=ceil((GPSstart+dur-POINT_SOURCE_MEMORY.bufstart-POINT_SOURCE_MEMORY.bufdur)/POINT_SOURCE_MEMORY.bufspline);
    lost=POINT_SOURCE_MEMORY.bufdur + nspline*POINT_SOURCE_MEMORY.bufspline - POINT_SOURCE_MEMORY.bufforget;
    if lost > 0
      % crop the buffer
      indi=1 + lost*POINT_SOURCE_MEMORY.fsample;
      POINT_SOURCE_MEMORY.bufstart = POINT_SOURCE_MEMORY.bufstart+lost;
      POINT_SOURCE_MEMORY.bufdur   = POINT_SOURCE_MEMORY.bufdur-lost;
      POINT_SOURCE_MEMORY.buf1     = POINT_SOURCE_MEMORY.buf1(indi:end);
      POINT_SOURCE_MEMORY.buf2     = POINT_SOURCE_MEMORY.buf2(indi:end);
      %fprintf('Buffer cropped\n');
    end
  else
    % the buffer is useless - delete it
    POINT_SOURCE_MEMORY.bufstart = GPSstart;
    POINT_SOURCE_MEMORY.bufdur   = -POINT_SOURCE_MEMORY.bufspline; % indicates that there is also no rolled iff data in buf
    POINT_SOURCE_MEMORY.buf1     = [];
    POINT_SOURCE_MEMORY.buf2     = [];
    %fprintf('Buffer erased\n');
  end
  % define sine and cosine arrays
  N=2*POINT_SOURCE_MEMORY.bufspline*POINT_SOURCE_MEMORY.fsample;
  s = sin(pi*[0:N/2-1]'/N);
  c = cos(pi*[0:N/2-1]'/N);
  % keep filling the buffer until enough data is in the can
  while GPSstart + dur > POINT_SOURCE_MEMORY.bufstart + POINT_SOURCE_MEMORY.bufdur , 
    % start time of new chunk
    chunkstart  =POINT_SOURCE_MEMORY.bufstart+POINT_SOURCE_MEMORY.bufdur;
    calibsec    =chunkstart+POINT_SOURCE_MEMORY.bufspline;
    siderealtime=GPStoGreenwichMeanSiderealTime(calibsec);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% calibrate                                                                   %%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % calculate transfer functions appropriate for this and half of the 
      % next segment (need both of these response functions for splicing)
      if ( strncmp(POINT_SOURCE_MEMORY.alphaBetaFile1,   'none', length(POINT_SOURCE_MEMORY.alphaBetaFile1))   | ...
           strncmp(POINT_SOURCE_MEMORY.calCavGainFile1,  'none', length(POINT_SOURCE_MEMORY.calCavGainFile1))  | ...
           strncmp(POINT_SOURCE_MEMORY.calResponseFile1, 'none', length(POINT_SOURCE_MEMORY.calResponseFile1)) )
        % the data is already calibrated
        transfer1 = constructFreqSeries(ones(POINT_SOURCE_MEMORY.numFreqs,1), POINT_SOURCE_MEMORY.flow, POINT_SOURCE_MEMORY.deltaF);
      else
        [R1, responseOK1] = ...
          calculateResponse(POINT_SOURCE_MEMORY.t1,    POINT_SOURCE_MEMORY.f1,...
			    POINT_SOURCE_MEMORY.R01,   POINT_SOURCE_MEMORY.C01,...
			    POINT_SOURCE_MEMORY.alpha1,POINT_SOURCE_MEMORY.gamma1,...
			    calibsec,                  POINT_SOURCE_MEMORY.ASQchannel1);
        % if response function is bad, set flag and exit loop
        if (responseOK1 == false)
          badResponse = true; h1=[]; h2=[];
          fprintf('bad response for detector 1 in Point Source simulation\n');
          return; % exit function
        end  
        % convert to transfer functions (units: counts/strain) 
        transfer1 = convertResponse(POINT_SOURCE_MEMORY.f1,    R1, POINT_SOURCE_MEMORY.flow,...
				    POINT_SOURCE_MEMORY.deltaF,    POINT_SOURCE_MEMORY.numFreqs, 1, 0);
      end

      if ( strncmp(POINT_SOURCE_MEMORY.alphaBetaFile2,   'none', length(POINT_SOURCE_MEMORY.alphaBetaFile2))   | ...
           strncmp(POINT_SOURCE_MEMORY.calCavGainFile2,  'none', length(POINT_SOURCE_MEMORY.calCavGainFile2))  | ...
           strncmp(POINT_SOURCE_MEMORY.calResponseFile2, 'none', length(POINT_SOURCE_MEMORY.calResponseFile2)) )

        % the data is already calibrated
        transfer2 = constructFreqSeries(ones(POINT_SOURCE_MEMORY.numFreqs,1), POINT_SOURCE_MEMORY.flow, POINT_SOURCE_MEMORY.deltaF);
      else
        [R2, responseOK2] = ...
          calculateResponse(POINT_SOURCE_MEMORY.t2,    POINT_SOURCE_MEMORY.f2,...
			    POINT_SOURCE_MEMORY.R02,   POINT_SOURCE_MEMORY.C02,...
			    POINT_SOURCE_MEMORY.alpha2,POINT_SOURCE_MEMORY.gamma2,...
			    calibsec,                  POINT_SOURCE_MEMORY.ASQchannel2);
        % if response function is bad, set flag and exit loop
        if (responseOK2 == false)
          badResponse = true; h1=[]; h2=[];
          fprintf('bad response for detector 2 in Point Source simulation\n');
          return; % exit function
        end  
        % convert to transfer functions (units: counts/strain) 
        transfer2 = convertResponse(POINT_SOURCE_MEMORY.f2,    R2, POINT_SOURCE_MEMORY.flow,...
				    POINT_SOURCE_MEMORY.deltaF,    POINT_SOURCE_MEMORY.numFreqs, 1, 0);
       end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% end of calibration                                                          %%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get new data;
    fprintf('Simulating data between %d and %d\n',chunkstart,chunkstart+POINT_SOURCE_MEMORY.bufspline*2);
    [d1,d2] =  simulatePointSource(POINT_SOURCE_MEMORY.bufspline*2,...
                                   POINT_SOURCE_MEMORY.fsample,...
			           POINT_SOURCE_MEMORY.Hf,...
			           POINT_SOURCE_MEMORY.intLog,...
			           POINT_SOURCE_MEMORY.det1,...
			           POINT_SOURCE_MEMORY.det2,...
			           siderealtime,...
			           POINT_SOURCE_MEMORY.ra,...
			           POINT_SOURCE_MEMORY.decl,...
				   POINT_SOURCE_MEMORY.power,...
			           transfer1,transfer2,...
				   POINT_SOURCE_MEMORY.MakeIncoherent);
    if or(any(isnan(d1)),any(isnan(d2)))
      error('simulated point source data contains NaN');
    end
    % refill bufs with new data;
    if POINT_SOURCE_MEMORY.bufdur < 0
      %special case: buffer has been deleted
      indi=1;
      indf=N/2;
      POINT_SOURCE_MEMORY.buf1(indi:indf,1)= c.*d1(N/2+1:N);
      POINT_SOURCE_MEMORY.buf2(indi:indf,1)= c.*d2(N/2+1:N);
    else
      %regular case: still data in the buffer
      indi=1 + POINT_SOURCE_MEMORY.bufdur * POINT_SOURCE_MEMORY.fsample ;
      indf=indi + N/2 - 1;
      POINT_SOURCE_MEMORY.buf1(indi:indf,1)=POINT_SOURCE_MEMORY.buf1(indi:indf) + s.*d1(1:N/2);
      POINT_SOURCE_MEMORY.buf2(indi:indf,1)=POINT_SOURCE_MEMORY.buf2(indi:indf) + s.*d2(1:N/2);
      indi=indi+N/2;
      indf=indf+N/2;
      POINT_SOURCE_MEMORY.buf1(indi:indf,1)= c.*d1(N/2+1:N);
      POINT_SOURCE_MEMORY.buf2(indi:indf,1)= c.*d2(N/2+1:N);
    end;
    % update pointer
    POINT_SOURCE_MEMORY.bufdur=POINT_SOURCE_MEMORY.bufdur + POINT_SOURCE_MEMORY.bufspline;
  end;
end;
% requested data is now in the buffer, return it
indi=1 + (GPSstart-POINT_SOURCE_MEMORY.bufstart)*POINT_SOURCE_MEMORY.fsample;
indf=indi + dur*POINT_SOURCE_MEMORY.fsample - 1;
h1=POINT_SOURCE_MEMORY.buf1(indi:indf);
h2=POINT_SOURCE_MEMORY.buf2(indi:indf);

return;

