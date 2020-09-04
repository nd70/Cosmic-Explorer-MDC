function [adcdata,fft_values,f]=detectorNoise(adcdata,DetectorNoiseFile)
%
%  detectorNoise --- simulates time-domain data from a given noise curve
%
%  arguments: 
%             adcdata   - data structure in preproc
%             DetectorNoiseFile - file containing noise curve
%
%  Routine written by Michael Coughlin, Eric Thrane.
%  Contact coughlim@carleton.edu, ethrane@physics.umn.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% convert sample rate and duration to number of samples and sample period
duration=length(adcdata.data)*adcdata.deltaT;
deltaT1=adcdata.deltaT;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

deltaT = deltaT1;
deltaF = 1/duration;
N=1/(deltaT*deltaF);

if ( mod(N,2)== 0 )
  numFreqs = N/2 - 1;
else
  numFreqs = (N-1)/2;
end

f = deltaF*[1:1:numFreqs]';
flow = deltaF;

NFFT = 2^nextpow2(N); % Next power of 2 from length of y

transfer=load(DetectorNoiseFile);
amp_values=transfer(:,2);
f_transfer1=transfer(:,1);
Pf1 = interp1(f_transfer1,amp_values, f, 'spline');

norm1 = sqrt(N/(2*deltaT)) * sqrt(Pf1);

re1 = norm1*sqrt(1/2).* randn(numFreqs, 1);
im1 = norm1*sqrt(1/2).* randn(numFreqs, 1);
z1  = re1 + 1i*im1;

% freq domain solution for htilde1, htilde2 in terms of z1, z2
htilde1 = z1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: what's below is taken from simulateSB.m isotropic code

% convolve data with instrument transfer function
otilde1 = htilde1*1;

% set DC and Nyquist = 0, then add negative freq parts in proper order
if ( mod(N,2)==0 )
  % note that most negative frequency is -f_Nyquist when N=even 
  otilde1 = [ 0; otilde1; 0; flipud(conj(otilde1)) ];
else
  % no Nyquist frequency when N=odd
  otilde1 = [ 0; otilde1; flipud(conj(otilde1)) ];
end;

% fourier transform back to time domain
o1_data = ifft(otilde1);

% take real part (imag part = 0 to round-off)
o1_data = real(o1_data);

adcdata.data=o1_data;

% psd estimation parameters
%Nt = length(d1);
%nfft = Nt/160;
%window = hann(nfft);
%noverlap = nfft/2;

% estimate power spectra from time-series data
%[p1,f1]=pwelch(d1,window,noverlap,nfft,fsample);
%[p2,f2]=pwelch(d2,window,noverlap,nfft,fsample);

Y = fft(adcdata.data,NFFT)/N;
fft_values=2*(abs(Y(1:NFFT/2-1)).^2)*duration;


return
