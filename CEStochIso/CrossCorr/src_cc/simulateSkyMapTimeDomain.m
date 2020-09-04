function [d1, d2] = simulateSkyMapTimeDomain(duration, tlow, ...
                                             fsample1, fsample2, ...
                                             nResample1, nResample2, betaParam1, betaParam2, ...
                                             Hf, intLog, ...
                                             siderealtime, ...
                                             glm, g1lm, g2lm, ...
                                             det1, det2, ...
                                             isSpH, coord1, coord2, map, ...
                                             transfer1, transfer2, ...
                                             seed)
%
%  simulateSkyMapTimeDomain --- simulates time-domain data for a spatially-distributed unpolarized SGWB
%
%  arguments: 
%             duration   - duration of the time series in seconds
%             tlow       - initial GPS time
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
%             siderealtime - sidereal time in hours (0..24h)
%             glm, g1lm, g2lm - spherical harmonic components of the overlap reduction functions (for detectors 12, 11, 22)
%             det1,det2  - detector structures containing position r and tensor d
%             isSpH      - Type of map: true: map contains complex spherical harmonics; false: map is pixel map
%             coord1     - either vector of l or right ascension in hours
%                          takes vector for multiple point sources, must have same length as coord2
%             coord2     - either vector of m or declination in degrees of source in the sky
%                          takes vector for multiple point sources, must have same length as coord1
%             map        - map data; either complex spherical haromnics, or value at pixel
%             transfer1,2- transfer function of detector 1/2
%
%  output:    d1,d2      - time-series data seen by detector 1,2
%
%  Routine written by Stefan Ballmer, Joe Romano.
%  Contact sballmer@ligo.mit.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set seed, if passed as an argument
try
  seed;
  randn('state',seed);
catch;
  % do nothing
end;

% convert sample rate and duration to number of samples and sample period
N1=floor(duration*fsample1);
N2=floor(duration*fsample2);

deltaT1=1/fsample1;
deltaT2=1/fsample2;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% discrete positive frequencies (not including DC component and Nyquist)
if N1 >= N2 
  N = N1;
  deltaT = deltaT1;
else
  N = N2;
  deltaT = deltaT2;
end

if ( mod(N,2)== 0 )
  numFreqs = N/2 - 1;
else
  numFreqs = (N-1)/2;
end

deltaF = 1/duration;
f = deltaF*[1:1:numFreqs]';
flow = deltaF;

% construct expected signal power/CP covariance matrix
if isSpH == true
  % convert Hf into a frequency series appropriate for SpH2CrossPower and Map2CrossPower
  H=checkArgumentFreqSeries(Hf,glm.flow,glm.deltaF,glm.numFreqs,intLog);

  % map is complex plms
  lvec = coord1;
  mvec = coord2;
  plm = map;

  % calculate expected signal power in an individual detector and 
  % the expected cross-power for the given plm's
  % (covariance matrix of htilde1, htilde2 is C = [ a c ; c b ])
  a = SpH2CrossPower(g1lm, plm, siderealtime, ...
                     numFreqs, flow, deltaF, H.data);
  a = real(a);
  b = SpH2CrossPower(g2lm, plm, siderealtime, ...
                     numFreqs, flow, deltaF, H.data);
  b = real(b);
  c = SpH2CrossPower(glm, plm, siderealtime, ...
                     numFreqs, flow, deltaF, H.data);
else
  % convert Hf into a frequency series appropriate for SpH2CrossPower and Map2CrossPower
  H=checkArgumentFreqSeries(Hf,flow,deltaF,numFreqs,intLog);

  % pixel map
  ra = coord1;
  decl = coord2;
  RADECL = [ra decl];
  dOmega = 1;

  % calculate expected signal power in an individual detector and 
  % the expected cross-power for the given map
  % (covariance matrix of htilde1, htilde2 is C = [ a c ; c b ])
  a = Map2CrossPower(det1, det1, map, dOmega, RADECL, siderealtime, ...
                     numFreqs, flow, deltaF, H.data);
  a = real(a);
  b = Map2CrossPower(det2, det2, map, dOmega, RADECL, siderealtime, ...
                     numFreqs, flow, deltaF, H.data);
  b = real(b);
  c = Map2CrossPower(det1, det2, map, dOmega, RADECL, siderealtime, ...
                     numFreqs, flow, deltaF, H.data);

end

% print out warning message if a, b, c correspond to an unphysical signal
detC = a.*b-c.*conj(c);
if min(detC) <= 0
  warning('Trying to simulate an unphysical signal')
end

% calculate eigenvalues of covariance matrix C = [ a c ; c b ]
r = sqrt((a-b).^2 + 4*c.*conj(c));
R = sqrt((r-a+b).^2 + 4*c.*conj(c));
lambda1 = 0.5*(a+b+r);
lambda2 = 0.5*(a+b-r);

% construct real and imaginary parts of uncorrelated random variables, 
% with random phases
norm = sqrt(N/(2*deltaT));

re1 = norm*sqrt(1/2).* randn(numFreqs, 1);
im1 = norm*sqrt(1/2).* randn(numFreqs, 1);
z1  = re1 + 1i*im1;

re2 = norm*sqrt(1/2).* randn(numFreqs, 1);
im2 = norm*sqrt(1/2).* randn(numFreqs, 1);
z2  = re2 + 1i*im2;

% freq domain solution for htilde1, htilde2 in terms of z1, z2
htilde1 = (2*conj(c).*sqrt(lambda1).*z1 + (r-a+b).*sqrt(lambda2).*z2)./R;
htilde2 = ((r-a+b).*sqrt(lambda1).*z1 - 2*c.*sqrt(lambda2).*z2)./R;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: what's below is taken from simulateSB.m isotropic code

% interpolate detector1 transfer functions to appropriate frequencies
f_transfer1 = transfer1.flow + transfer1.deltaF*[0:length(transfer1.data)-1]';
f_transfer2 = transfer2.flow + transfer2.deltaF*[0:length(transfer2.data)-1]';
transfer1_i = interp1(f_transfer1, transfer1.data, f, 'spline');
transfer2_i = interp1(f_transfer2, transfer2.data, f, 'spline');

% convolve data with instrument transfer function
otilde1 = htilde1.*transfer1_i;
otilde2 = htilde2.*transfer2_i;

% set DC and Nyquist = 0, then add negative freq parts in proper order
if ( mod(N,2)==0 )
  % note that most negative frequency is -f_Nyquist when N=even 
  otilde1 = [ 0; otilde1; 0; flipud(conj(otilde1)) ];
  otilde2 = [ 0; otilde2; 0; flipud(conj(otilde2)) ];
else
  % no Nyquist frequency when N=odd
  otilde1 = [ 0; otilde1; flipud(conj(otilde1)) ];
  otilde2 = [ 0; otilde2; flipud(conj(otilde2)) ];
end;

% fourier transform back to time domain
o1_data = ifft(otilde1);
o2_data = ifft(otilde2);

% take real part (imag part = 0 to round-off)
o1_data = real(o1_data);
o2_data = real(o2_data);

if ( N1 ~= N )
  o1_data = resample(o1_data, N1, N, nResample1, betaParam1);
end

if ( N2 ~= N )
  o2_data = resample(o2_data, N2, N, nResample2, betaParam2);
end

d1 = o1_data;
d2 = o2_data;

return
