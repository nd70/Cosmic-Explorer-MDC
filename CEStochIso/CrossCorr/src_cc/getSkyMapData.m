function [h1, h2,badResponse] = getSkyMapData(GPSstart,dur)
%  function [h1, h2] = getSkyMapData(GPSstart,dur)
%
%  getSkyMapData --- returns dur seconds of simulated sky map data starting at GPSstart
%                    The function relies on the GLOBAL variable SKY_MAP_MEMORY,
%                    which is a struct and has to be initialized properly before the 1st call.
%
%  arguments: GPSstart   - GPS start time of the data
%             dur        - duration of the data in seconds
%
%  SKY_MAP_MEMORY:
%             injectTimeDomain - if true, inject in time-domain; if false, inject in frequency domain
%             fsample1   - sample frequency of the time series 1
%             fsample2   - sample frequency of the time series 2
%             nResample1 - order of matlab resample routine for time series 1
%             nResample2 - order of matlab resample routine for time series 2
%             betaParam1 - beta parameter for matlab resample routine for time series 1
%             betaParam2 - beta parameter for matlab resample routine for time series 2
%             Hf         - total power spectrum (one-sided) in both polarizations, i.e.
%                          the actual power spectrum for each polarization is Hf/2
%                          preferred input: Nx2 array (filename and freq. series
%                          work too ... in principle)
%             intLog     - boolean whether to interpolate Hf logarithmicly
%             Lmax       - Lmax for simulated sky map if isSpH is true (may differ from the Lmax used for analysis)
%             ifo1, ifo2 - strings identifying the two detectors (e.g., H1, L1)
%             glm, g1lm, g2lm - spherical harmonic components of the overlap reduction functions (for detectors 12, 11, 22)
%             det1,det2  - detector structures containing position r and tensor d
%             isSpH      - Type of map: true: map contains complex spherical harmonics; false: map is pixel map
%             coord1     - either vector of l or right ascension in hours
%                          takes vector for multiple point sources, must have same length as coord2
%             coord2     - either vector of m or declination in degrees of source in the sky
%                          takes vector for multiple point sources, must have same length as coord1
%	      map        - map data; either complex spherical haromnics, or value at pixel
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
%  Routine written by Stefan Ballmer, Joe Romano.
%  Contact sballmer@ligo.mit.edu
%
%  $Id: getSkyMapData.m,v 1.13 2008-09-25 17:30:52 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global SKY_MAP_MEMORY;
badResponse = false;

% is GPSstart too far back?
if GPSstart < SKY_MAP_MEMORY.bufstart
  error('cannot go that far into the past - increase bufforget');
end

% Do I need to get new data?
if GPSstart + dur > SKY_MAP_MEMORY.bufstart + SKY_MAP_MEMORY.bufdur
  % does the buffer overlap with the new past?
  if GPSstart+dur- SKY_MAP_MEMORY.bufforget < ...
                   SKY_MAP_MEMORY.bufstart + SKY_MAP_MEMORY.bufdur + SKY_MAP_MEMORY.bufspline
    % does the buffer get too long?
    nspline=ceil((GPSstart+dur-SKY_MAP_MEMORY.bufstart-SKY_MAP_MEMORY.bufdur)/SKY_MAP_MEMORY.bufspline);
    lost=SKY_MAP_MEMORY.bufdur + nspline*SKY_MAP_MEMORY.bufspline - SKY_MAP_MEMORY.bufforget;
    if lost > 0
      % crop the buffer
      indi1=1 + lost*SKY_MAP_MEMORY.fsample1;
      indi2=1 + lost*SKY_MAP_MEMORY.fsample2;
      SKY_MAP_MEMORY.bufstart = SKY_MAP_MEMORY.bufstart+lost;
      SKY_MAP_MEMORY.bufdur   = SKY_MAP_MEMORY.bufdur-lost;
      SKY_MAP_MEMORY.buf1     = SKY_MAP_MEMORY.buf1(indi1:end);
      SKY_MAP_MEMORY.buf2     = SKY_MAP_MEMORY.buf2(indi2:end);
      %fprintf('Buffer cropped\n');
    end
  else
    % the buffer is useless - delete it
    SKY_MAP_MEMORY.bufstart = GPSstart;
    SKY_MAP_MEMORY.bufdur   = -SKY_MAP_MEMORY.bufspline; % indicates that there is also no rolled iff data in buf
    SKY_MAP_MEMORY.buf1     = [];
    SKY_MAP_MEMORY.buf2     = [];
    %fprintf('Buffer erased\n');
  end
  % define sine and cosine arrays
  N1=2*SKY_MAP_MEMORY.bufspline*SKY_MAP_MEMORY.fsample1;
  N2=2*SKY_MAP_MEMORY.bufspline*SKY_MAP_MEMORY.fsample2;
  s1 = sin(pi*[0:N1/2-1]'/N1);
  s2 = sin(pi*[0:N2/2-1]'/N2);
  c1 = cos(pi*[0:N1/2-1]'/N1);
  c2 = cos(pi*[0:N2/2-1]'/N2);
  % keep filling the buffer until enough data is in the can
  while GPSstart + dur > SKY_MAP_MEMORY.bufstart + SKY_MAP_MEMORY.bufdur , 
    % start time of new chunk
    chunkstart  =SKY_MAP_MEMORY.bufstart+SKY_MAP_MEMORY.bufdur;
    calibsec    =chunkstart+SKY_MAP_MEMORY.bufspline;
    siderealtime=GPStoGreenwichMeanSiderealTime(calibsec);
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% calibrate                                                                   %%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % calculate transfer functions appropriate for this and half of the 
      % next segment (need both of these response functions for splicing)
      if ( strncmp(SKY_MAP_MEMORY.alphaBetaFile1,   'none', length(SKY_MAP_MEMORY.alphaBetaFile1))   | ...
           strncmp(SKY_MAP_MEMORY.calCavGainFile1,  'none', length(SKY_MAP_MEMORY.calCavGainFile1))  | ...
           strncmp(SKY_MAP_MEMORY.calResponseFile1, 'none', length(SKY_MAP_MEMORY.calResponseFile1)) )
        % the data is already calibrated
        transfer1 = constructFreqSeries(ones(SKY_MAP_MEMORY.numFreqs,1), SKY_MAP_MEMORY.flow, SKY_MAP_MEMORY.deltaF);
      else
        [R1, responseOK1] = ...
          calculateResponse(SKY_MAP_MEMORY.t1,    SKY_MAP_MEMORY.f1,...
			    SKY_MAP_MEMORY.R01,   SKY_MAP_MEMORY.C01,...
			    SKY_MAP_MEMORY.alpha1,SKY_MAP_MEMORY.gamma1,...
			    calibsec,                  SKY_MAP_MEMORY.ASQchannel1);
        % if response function is bad, set flag and exit loop
        if (responseOK1 == false)
          badResponse = true; h1=[]; h2=[];
          fprintf('bad response for detector 1 in Point Source simulation\n');
          return; % exit function
        end  
        % convert to transfer functions (units: counts/strain) 
        transfer1 = convertResponse(SKY_MAP_MEMORY.f1,    R1, SKY_MAP_MEMORY.flow,...
				    SKY_MAP_MEMORY.deltaF,    SKY_MAP_MEMORY.numFreqs, 1, 0);
      end

      if ( strncmp(SKY_MAP_MEMORY.alphaBetaFile2,   'none', length(SKY_MAP_MEMORY.alphaBetaFile2))   | ...
           strncmp(SKY_MAP_MEMORY.calCavGainFile2,  'none', length(SKY_MAP_MEMORY.calCavGainFile2))  | ...
           strncmp(SKY_MAP_MEMORY.calResponseFile2, 'none', length(SKY_MAP_MEMORY.calResponseFile2)) )

        % the data is already calibrated
        transfer2 = constructFreqSeries(ones(SKY_MAP_MEMORY.numFreqs,1), SKY_MAP_MEMORY.flow, SKY_MAP_MEMORY.deltaF);
      else
        [R2, responseOK2] = ...
          calculateResponse(SKY_MAP_MEMORY.t2,    SKY_MAP_MEMORY.f2,...
			    SKY_MAP_MEMORY.R02,   SKY_MAP_MEMORY.C02,...
			    SKY_MAP_MEMORY.alpha2,SKY_MAP_MEMORY.gamma2,...
			    calibsec,             SKY_MAP_MEMORY.ASQchannel2);
        % if response function is bad, set flag and exit loop
        if (responseOK2 == false)
          badResponse = true; h1=[]; h2=[];
          fprintf('bad response for detector 2 in Point Source simulation\n');
          return; % exit function
        end  
        % convert to transfer functions (units: counts/strain) 
        transfer2 = convertResponse(SKY_MAP_MEMORY.f2,    R2, SKY_MAP_MEMORY.flow,...
				    SKY_MAP_MEMORY.deltaF,    SKY_MAP_MEMORY.numFreqs, 1, 0);
       end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%% end of calibration                                                          %%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get new data;
    fprintf('Simulating sky map data between %d and %d\n',chunkstart,chunkstart+SKY_MAP_MEMORY.bufspline*2);
    [d1,d2] =  simulateSkyMapTimeDomain(SKY_MAP_MEMORY.bufspline*2,...
                                        GPSstart,...
                                        SKY_MAP_MEMORY.fsample1,...
                                        SKY_MAP_MEMORY.fsample2,...
                                        SKY_MAP_MEMORY.nResample1,...
                                        SKY_MAP_MEMORY.nResample2,...
                                        SKY_MAP_MEMORY.betaParam1,...
                                        SKY_MAP_MEMORY.betaParam2,...
			                SKY_MAP_MEMORY.Hf,...
                                        SKY_MAP_MEMORY.intLog,...
			                siderealtime,...
                                        SKY_MAP_MEMORY.glm,...
                                        SKY_MAP_MEMORY.g1lm,...
                                        SKY_MAP_MEMORY.g2lm,...
                                        SKY_MAP_MEMORY.det1,...
                                        SKY_MAP_MEMORY.det2,...
                                        SKY_MAP_MEMORY.isSpH,...
			                SKY_MAP_MEMORY.coord1,...
			                SKY_MAP_MEMORY.coord2,...
		                        SKY_MAP_MEMORY.map,...
                                        transfer1,transfer2);
    if or(any(isnan(d1)),any(isnan(d2)))
      error('simulated sky map data contains NaN');
    end
    % refill bufs with new data;
    if SKY_MAP_MEMORY.bufdur < 0
      %special case: buffer has been deleted
      indi1=1;
      indi2=1;
      indf1=N1/2;
      indf2=N2/2;
      SKY_MAP_MEMORY.buf1(indi1:indf1,1)= c1.*d1(N1/2+1:N1);
      SKY_MAP_MEMORY.buf2(indi2:indf2,1)= c2.*d2(N2/2+1:N2);
    else
      %regular case: still data in the buffer
      indi1=1 + SKY_MAP_MEMORY.bufdur * SKY_MAP_MEMORY.fsample1;
      indi2=1 + SKY_MAP_MEMORY.bufdur * SKY_MAP_MEMORY.fsample2;
      indf1=indi1 + N1/2 - 1;
      indf2=indi2 + N2/2 - 1;
      SKY_MAP_MEMORY.buf1(indi1:indf1,1)=SKY_MAP_MEMORY.buf1(indi1:indf1) + s1.*d1(1:N1/2);
      SKY_MAP_MEMORY.buf2(indi2:indf2,1)=SKY_MAP_MEMORY.buf2(indi2:indf2) + s2.*d2(1:N2/2);
      indi1=indi1+N1/2;
      indi2=indi2+N2/2;
      indf1=indf1+N1/2;
      indf2=indf2+N2/2;
      SKY_MAP_MEMORY.buf1(indi1:indf1,1)= c1.*d1(N1/2+1:N1);
      SKY_MAP_MEMORY.buf2(indi2:indf2,1)= c2.*d2(N2/2+1:N2);
    end;
    % update pointer
    SKY_MAP_MEMORY.bufdur=SKY_MAP_MEMORY.bufdur + SKY_MAP_MEMORY.bufspline;
  end;
end;
% requested data is now in the buffer, return it
indi1=1 + (GPSstart-SKY_MAP_MEMORY.bufstart)*SKY_MAP_MEMORY.fsample1;
indi2=1 + (GPSstart-SKY_MAP_MEMORY.bufstart)*SKY_MAP_MEMORY.fsample2;
indf1=indi1 + dur*SKY_MAP_MEMORY.fsample1 - 1;
indf2=indi2 + dur*SKY_MAP_MEMORY.fsample2 - 1;
h1=SKY_MAP_MEMORY.buf1(indi1:indf1);
h2=SKY_MAP_MEMORY.buf2(indi2:indf2);

return;

