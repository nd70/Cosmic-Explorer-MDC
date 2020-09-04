function [h1, h2] = getPointSourceData_IM(injType, GPSstart, dur, ...
  pp, fsample, Hf, det1, det2, ra, decl, power_amplifier, MakeIncoherent, ...
  intLog, params)
%
%  getPointSourceData --- returns dur seconds of simulated point source data 
%  starting at GPSstart
%
%  arguments: GPSstart   - injType 
%                          'PSD' for power spectral injection
%                          'time-series' for time-series injection
%             GPSstart   - GPSt start time of the signal
%             dur        - duration of the signal in seconds
%                          if dur==-1, dur is resent by file length
%             fsample    - sample frequency of the time series
%             Hf         - i) total power spectrum (one-sided) in both
%                            polarizations, i.e.
%                     the actual power spectrum for each polarization is Hf/2
%                          input: Nx2 array (freq. series and PSD)
%                      ii) in case of 'timeseries' injection, it will be Nx2
%                          array (time stamps and h(t))
%
%          det1,det2  - detector structures containing position r and tensor d
%             ra         - right ascension in hours of source in the sky
%     takes vector for multiple point sources, must have same length as decl
%             decl       - right ascension in degrees of source in the sky
%         takes vector for multiple point sources, must have same length as ra
%             pp.flow       \      
%             pp.deltaF      > frequency information used for the main analysis
%             pp.fhigh      /   
%      MakeIncoherent optional parameter to destroy coherence of point source
%                            0:  coherent point source
%                        1:  incoherent, but scaled with DC antenna acceptance
%        2:  stationary noise, PowerSpec corresponds to the noise seen in both 
%
%  output:    h1,h2   - time series seen by detector 1,2, calibrated in counts
%
%  Routine adopted from getPointSourcedata.m written by Stefan Ballmer.
%  Contact shivaraj@physics.umn.edu, ethrane
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numFreqs = floor((pp.fhigh-pp.flow)/pp.deltaF)+1;
siderealtime = GPStoGreenwichMeanSiderealTime(GPSstart);

% The data is already calibrated and hence the transfer function is one 
transfer1 = constructFreqSeries(ones(numFreqs,1), pp.flow, pp.deltaF);
transfer2 = constructFreqSeries(ones(numFreqs,1), pp.flow, pp.deltaF);

if(strcmp(injType,'PSD'))
  % Simulating the data;
  fprintf('Simulating data between %d and %d\n',GPSstart,GPSstart+dur);

  [d1,d2] =  simulatePointSource(dur,fsample,Hf,intLog,det1,det2, ...
              siderealtime, ra,decl,power_amplifier,transfer1,...
                                  transfer2,MakeIncoherent);

  if or(any(isnan(d1)),any(isnan(d2)))
     error('Simulated point source data contains NaN.\n');
  end

  h1=d1;
  h2=d2;
end

if(strcmp(injType,'time_series'))
  if(size(Hf,2)==3)
    hp = Hf(:,2);
    hx = Hf(:,3); 
  else
    hp = Hf(:,2);
    hx = Hf(:,2);
  end

  if dur==-1
    dur = Hf(end,1);
  end

  tt = 1/fsample:1/fsample:dur;
  hp = interp1(Hf(:,1)- Hf(1,1),hp,tt,'spline')';
  hx = interp1(Hf(:,1)- Hf(1,1),hx,tt,'spline')';

  h_length=fsample*dur;
  h1=zeros(h_length,1);
  h2=zeros(h_length,1);

  for kk=1:length(ra)
    g=orfIntegrandSymbolic(det1,det2,siderealtime,ra(kk),decl(kk));
    try
      if params.fixAntennaFactors
        g=fixAntennaFactors(g);
      end
    catch
    end
    switch MakeIncoherent
       case 0  % 0:  coherent point source
%fprintf('############ THIS IS A TEST ############\n');
%fprintf('############ THIS IS A TEST ############\n');
%fprintf('############ THIS IS A TEST ############\n');
%g.F1p=1;
%g.F2p=1;
           h1_ra_dec= (hp*g.F1p + hx*g.F1x) * sqrt(power_amplifier(kk));
           h2_ra_dec= (hp*g.F2p + hx*g.F2x) * sqrt(power_amplifier(kk));
       case 1  % 1:  incoherent, but scaled with DC antenna acceptance
           h1_ra_dec= hp*sqrt((g.F1p^2 + g.F1x^2) * power_amplifier(kk));
           h2_ra_dec= hx*sqrt((g.F2p^2 + g.F2x^2) * power_amplifier(kk));
       case 2  % 2:  stationary noise, 
               % PowerSpec corresponds to the noise seen in both
           h1_ra_dec=hp*sqrt(2 * power_amplifier(kk));
           h2_ra_dec=hx*sqrt(2 * power_amplifier(kk));
       otherwise
           error('Unknown value of MakeIncoherent in simulatePointSource.\n');
    end

    % make h1 have the same array size as h2 by
    % fourier transforming it and ifft-ing it back
    h1_ra_dec_shift=shiftTimeSeries(h1_ra_dec,fsample,0);

    % fourier transform h2 to create h2_tilde
    % multiply h0_tilde by exp(2pi*i*f*tau) to transform h0(t) to h0(t+tau)
    h2_ra_dec_shift=shiftTimeSeries(h2_ra_dec,fsample,g.tau);
	
    %plot time shifted series and un-timeshifted series
    h1=h1+h1_ra_dec_shift;
    h2=h2+h2_ra_dec_shift;

  end %end loop over RA and dec
end

return;

