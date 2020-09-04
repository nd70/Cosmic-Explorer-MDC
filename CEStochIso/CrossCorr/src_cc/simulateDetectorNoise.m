function [d1, d2] = simulateDetectorNoise(duration, tlow, ...
                                          fsample1, fsample2, ...
                                          nResample1, nResample2, betaParam1, betaParam2, ...
                                          P1, P2, intLog, ...
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
%             P1,2       - power spectra (one-sided, strain^2/Hz) in detectors 1,2
%                          preferred input: Nx2 array (filename and freq. series
%                          work too ... in principle)
%             intLog     - boolean whether to interpolate P1,2 logarithmicly
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

% convert P1,2 to a frequency series structure
Pf1=checkArgumentFreqSeries(P1,flow,deltaF,numFreqs,intLog);
Pf2=checkArgumentFreqSeries(P2,flow,deltaF,numFreqs,intLog);

% construct real and imaginary parts of uncorrelated random variables, 
% with random phases, weighted by respective amplitude spectral densities
norm1 = sqrt(N/(2*deltaT)) * sqrt(Pf1.data);
norm2 = sqrt(N/(2*deltaT)) * sqrt(Pf2.data);

re1 = norm1*sqrt(1/2).* randn(numFreqs, 1);
im1 = norm1*sqrt(1/2).* randn(numFreqs, 1);
z1  = re1 + 1i*im1;

re2 = norm2*sqrt(1/2).* randn(numFreqs, 1);
im2 = norm2*sqrt(1/2).* randn(numFreqs, 1);
z2  = re2 + 1i*im2;

% freq domain solution for htilde1, htilde2 in terms of z1, z2
htilde1 = z1;
htilde2 = z2;

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
