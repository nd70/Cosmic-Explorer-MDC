function [h1, h2] = simulatePointSource(dur,fsample,Hf,intLog,det1,det2,time,ra,decl,power,transfer1,transfer2,MakeIncoherent)
%  function [h1, h2] = simulatePointSource(dur,fsample,Hf,intLog,det1,det2,time,ra,decl,power,transfer1,transfer2,MakeIncoherent)
%
%  simulatePointSource --- simulates a stochastic unpolarized point source in the sky
%
%  arguments: dur        - duration of the time series in seconds
%             fsample    - sample frequency of the time series
%             Hf         - total power spectrum (one-sided) in both polarizations, i.e.
%                          the actual power spectrum for each polarization is Hf/2
%                          preferred input: Nx2 array (filename and freq. series
%                          work too ... in principle)
%             intLog     - boolean whether to interpolate Hf logarithmicly
%             det1,det2  - detector structures containing position r and tensor d
%             time       - sidereal time in hours (0..24h)
%             ra         - right ascension in hours of source in the sky
%                          takes vector for multiple point sources, must have same length as decl
%             decl       - right ascension in degrees of source in the sky
%                          takes vector for multiple point sources, must have same length as ra
%	      power      - power emanating from corresponding point source (aka coefficient of Hf for a given source)
%			   takes vector for multiple point sources, must have same length as ra and decl
%             transfer1,2- Transfer function of detector 1/2
%             MakeIncoherent optional parameter to destroy coherence of point source
%                            0:  coherent point source
%                            1:  incoherent, but scaled with DC antenna acceptance
%                            2:  stationary noise, PowerSpec corresponds to the noise seen in both 
%
%  output:    h1,h2      - time series seen by detector 1,2
%
%  Routine written by Stefan Ballmer.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
  MakeIncoherent;
catch
  MakeIncoherent = 0;
end;

doCal=true;

try
  transfer1;
catch
  doCal = false;
end;

try
  transfer2;
catch
  doCal = false;
end;

N=dur*fsample;
h1=zeros(N,1);
h2=zeros(N,1);

for kk=1:length(ra)
  g=orfIntegrandSymbolic(det1,det2,time,ra(kk),decl(kk));
  % tau = tarrival1 - tarrival2, i.e. tau > 0 means det 2 is closer to source
  if doCal
    [hp,hx]=randomTS(dur,fsample,Hf,intLog,2,[0,g.tau],[],transfer1,transfer2);
  else
    [hp,hx]=randomTS(dur,fsample,Hf,intLog,2,[0,g.tau]);
  end;

  switch MakeIncoherent
    case 0  % 0:  coherent point source
      h1=h1+ (hp(:,1)*g.F1p + hx(:,1)*g.F1x) * sqrt(power(kk));
      h2=h2+ (hp(:,2)*g.F2p + hx(:,2)*g.F2x) * sqrt(power(kk));
    case 1  % 1:  incoherent, but scaled with DC antenna acceptance
      h1=h1+ hp(:,1)*sqrt((g.F1p^2 + g.F1x^2) * power(kk));
      h2=h2+ hx(:,2)*sqrt((g.F2p^2 + g.F2x^2) * power(kk));
    case 2  % 2:  stationary noise, PowerSpec corresponds to the noise seen in both
      h1=h1+hp(:,1)*sqrt(2 * power(kk));
      h2=h2+hx(:,2)*sqrt(2 * power(kk));
    otherwise
      error('Unknown value of MakeIncoherent in simulatePointSource');
      return;
  end;
end;
return

