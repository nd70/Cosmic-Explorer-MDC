function [o1, o2] = simulateSB(tlow, deltaT1, deltaT2, N1, N2,...
                               signalType, detector1, detector2, ...
			       transfer1, transfer2, nResample1, betaParam1,...
                               nResample2, betaParam2, fbase1, fbase2, seed)
%
% simulateSB --- simulates a stochastic background signal
%
% simulateSB(tlow, deltaT1, deltaT2, N1, N2, signalType, detector1, detector2,
% transfer1, transfer2, nResample1, betaParam1, nResample2, betaParam2, fbase1,
% fbase2, seed) generates N1 and N2 samples (for detectors1,2, respectively) of
% two simulated unpolarized, isotropic stochastic background signals convolved 
% with the instrument transfer functions (units of o1,2: counts)
%
%  signalType:
%  - 'white' -> SB signal with S_gw(f) = (3 H_100^2/10 pi^2) 
%  - 'const' -> SB signal with h_100^2 omegaGW = 1
% 
% NOTE: the white case corresponds to h_100^2 omegaGW(f)/f^3 = 1/Hz^3
% or h_100^2 omegaGW(fRef) = (fRef/Hz)^3.
% 
% This routine works for any GW detector pair, namely, ifo-ifo, ifo-bar,
% bar-ifo, or bar-bar.
%
% It still requires testing for the case of nonzero tlow.
%
% transfer1,2 = instrument transfer functions (units: counts/strain)
%
% Routine written by Joseph D. Romano, Sukanta Bose, John T. Whelan
%  Contact Joseph.Romano@astro.cf.ac.uk%
%  $Id: simulateSB.m,v 1.8 2008-10-11 02:48:42 jromano Exp $
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set seed, if passed as an argument
try 
  seed;
  randn('state',seed);
catch
  % do nothing
end

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

deltaF = 1/(N*deltaT);
f = deltaF*[1:1:numFreqs]';

% normalization depends on signal type
if isnumeric(signalType) % signalType = powerIndex
  H100 = HubbleConstant;
  norm = sqrt(N/(2*deltaT)) * sqrt((3*H100^2)/(10*pi^2));
  % frequency dependence is f^-((powerIndex-3)/2)
  norm = norm*f.^((signalType-3)/2);
else 
  if(strncmp(signalType, 'white', length(signalType)))
    % S_gw(f) = (3 H_100^2/10 pi^2) 
    H100 = HubbleConstant;
    norm = sqrt(N/(2*deltaT)) * sqrt((3*H100^2)/(10*pi^2));
  elseif ( strncmp(signalType, 'const', length(signalType)) )
    % h_100^2 omegaGW(f) = 1 
    H100 = HubbleConstant;
    norm = sqrt(N/(2*deltaT)) * sqrt((3*H100^2)/(10*pi^2));
    % frequency dependence is f^-(3/2)
    norm = norm*f.^(-3/2);
  end
end

% construct real and imaginary parts, with random phases
re1 = (norm/sqrt(2)).* randn(numFreqs, 1);
im1 = (norm/sqrt(2)).* randn(numFreqs, 1);
z1  = re1 + 1i*im1;

re2 = (norm/sqrt(2)).* randn(numFreqs, 1);
im2 = (norm/sqrt(2)).* randn(numFreqs, 1);
z2  = re2 + 1i*im2;

% fourier components at site 2 are correlated with those at site 1
% via the overlap reduction functions          
% Generalized the way the strains are generated to account for possible 
%ifo-bar analysis
gamma12 = constructFreqSeries(overlapreductionfunction(f, detector1, detector2), ...
                            flow, deltaF);
gamma11 = constructFreqSeries(overlapreductionfunction(f, detector1, detector1), ...
                            flow, deltaF);
gamma22 = constructFreqSeries(overlapreductionfunction(f, detector2, detector2), ...
                            flow, deltaF);
sDenom = gamma11.data .* gamma22.data;
s = sqrt(1 - gamma12.data.^2 ./sDenom);

a = (1/sqrt(2)) * sqrt(1+s);
b = (1/sqrt(2)) * gamma12.data ./ (sqrt((1+s) .*sDenom));

htilde1 = sqrt(gamma11.data) .*(a.*z1 + b.*z2);
htilde2 = sqrt(gamma22.data) .*(b.*z1 + a.*z2);

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

if ( ~isnan(fbase1) )
  %heterodyne the waveform
  % t should really be tlow plus this, and phase should be included
  t = deltaT*[0:1:(N-1)]';
  o1_data = o1_data.*(cos(2*pi*fbase1.*t) - 1i*sin(2*pi*fbase1.*t));
end

if ( N1 ~= N )
  o1_data = resample(o1_data, N1, N, nResample1, betaParam1);
end

% fill structures for detector1
o1 = constructTimeSeries(o1_data, tlow, deltaT1);

if ( ~isnan(fbase2) )
  %heterodyne the waveform
  % t should really be tlow plus this, and phase should be included
  t = deltaT*[0:1:(N-1)]';
  o2_data = o2_data.*(cos(2*pi*fbase2.*t) - 1i*sin(2*pi*fbase2.*t));
end

if ( N2 ~= N )
  o2_data = resample(o2_data, N2, N, nResample2, betaParam2);
end

% fill structures for detector2
o2 = constructTimeSeries(o2_data, tlow, deltaT2);

return
