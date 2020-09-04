function [norm, theorVar, sensInt] = ...
  normalization(T, gamma, fRef, alpha, calP1, calP2, window1, window2, mask, pp)
%
%  normalization --- normalization for a stochastic optimal filter
%
%  normalization(T, gamma, fRef, alpha, calP1, calP2, window1, window2, mask,pp)
%  calculates the normalization factor for the optimal filter, the
%  theoretical variance of the optimally-filtered CC  statistic, and
%  the integrand of 1/(theoretical variance)--the so-called 
%  sensitivity integrand.  The normalization is chosen so that the 
%  expected value of the CC statistic is equal to omegaRef * T, where 
%  T is the time interval for a single stretch of analyzed data (e.g., 
%  60 s).  sensInt is a freq-series structure with associated metadata.
%
%  If neither time series was heterodyned (base-banded), then
%
%  norm      =  (10 \pi^2 / 3 H_100^2) * 2 *
%               (1 / mean(window1 * window2)) *
%               [2 * int_flow^fhigh deltaF gamma^2(f) (f/fRef)^(2*alpha) / 
%                f^6 P1(f) P2(f)]^-1
%
%  theorVar  =  T * (10 \pi^2 / 3 H_100^2)^2 *
%               mean(window1^2 * window2^2)/mean(window1 * window2)^2 *
%               [2 * int_flow^fhigh deltaF gamma^2(f) (f/fRef)^(2*alpha) /
%                f^6 P1(f) P2(f)]^-1
% 
%  sensInt   =  (1/T) * (3 H_100^2/ 10 \pi^2)^2 *
%               mean(window1 * window2)^2 / mean(window1^2 * window2^2) *
%               [2 * gamma^2(f) (f/fRef)^(2*alpha) /
%                f^6 P1(f) P2(f)]
%
%  If either time series was heterodyned (base-banded), then
%
%  norm      =  (10 \pi^2 / 3 H_100^2) * 2 *
%               (1 / mean(window1 * window2)) *
%               [ int_flow^fhigh deltaF gamma^2(f) (f/fRef)^(2*alpha) / 
%                f^6 P1(f) P2(f)]^-1
% 
%  theorVar  =  T * (10 \pi^2 / 3 H_100^2)^2 *
%               mean(window1^2 * window2^2)/mean(window1 * window2)^2 *
%               [ int_flow^fhigh deltaF gamma^2(f) (f/fRef)^(2*alpha) /
%                f^6 P1(f) P2(f)]^-1 / 2
% 
%  sensInt   =  (1/T) * (3 H_100^2/ 10 \pi^2)^2 *
%               mean(window1 * window2)^2 / mean(window1^2 * window2^2) *
%               [ gamma^2(f) (f/fRef)^(2*alpha) /
%                f^6 P1(f) P2(f)]
%
%  where theorVar is the theoretical variance of the real or imaginary
%  part of the complex cross-correlations statistic.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Extension for arbitrary signal power spectrum:
%
%  function [norm, theorVar, sensInt] = ...
%    normalization(T, gamma, Hf, [], calP1, calP2, window1, window2, mask)
%
%  i.e. instead of fref it takes a frequency series. This version is used
%  for the radiometer.
%
%  The normalization is basically the same as in the full-sky search, with
%  the following replacements:
%
%   -   3 H_100^2/ 10 \pi^2    ---->    1
%   -   (f/fRef)^alpha / f^3   ---->    Hf
%
%  This assumes that gamma=sum_A(F1_A * F2_A*exp(2*pi*i*f*tau)), i.e.
%  gamma(f=0) is no longer equal to 1 as in the full sky case.
%  The 5/(8*pi) in the gamma normalization was also dropped since it makes
%  no sense without integration over the sky.
%  Also note that H(f) is the total power summed over the 2 polarizations.
%
%  norm      =  2 *
%               (1 / mean(window1 * window2)) *
%               [2 * int_flow^fhigh deltaF abs(gamma)^2(f) Hf^2 / 
%                P1(f) P2(f)]^-1
%
%  theorVar  =  T *
%               mean(window1^2 * window2^2)/mean(window1 * window2)^2 *
%               [2 * int_flow^fhigh deltaF abs(gamma)^2(f) Hf^2 /
%                P1(f) P2(f)]^-1
% 
%  sensInt   =  (1/T) *
%               mean(window1 * window2)^2 / mean(window1^2 * window2^2) *
%               [2 * abs(gamma)^2(f) Hf^2 /
%                P1(f) P2(f)]
%
%  NOTE: mask zeroes out certain frequency bins, and pp is (optional) 
%  structure that contains the already calculated window parameters. If 
%  pp is specified, window1 and window2 should be set to -1.
% 
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org
%  Extendend by Stefan Ballmer, sballmer@ligo.mit.edu
% 
%  $Id: normalization.m,v 1.13 2008-09-12 15:54:32 ethrane Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract frequency series metadata
[data, flow, deltaF, symmetry] = extractFreqSeries(gamma);
numFreqs = length(data);
                                                                                
% if no mask, set to one
try
  mask;
catch
  mask = constructFreqSeries(ones(numFreqs,1), flow, deltaF);
end
                                                                                
% error checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that gamma, P1, P2, and mask have the same length
useHf = true;
try
 fRef.flow;
catch
  useHf = false;
end;

if ( (length(gamma.data) ~= numFreqs) | ...
     (length(calP1.data) ~= numFreqs) | ...
     (length(calP2.data) ~= numFreqs) | ...
     (length(mask.data)  ~= numFreqs) )
  error('size mismatch');
end

% check that gamma, P1, P2, and mask have the same flow
if ( (gamma.flow ~= flow) | ...
     (calP1.flow ~= flow) | ...
     (calP2.flow ~= flow) | ...
     (mask.flow  ~= flow) )
  error('flow mismatch');
end

% check that gamma, P1, P2, and mask have the same delatF
if ( (gamma.deltaF ~= deltaF) | ...
     (calP1.deltaF ~= deltaF) | ...
     (calP2.deltaF ~= deltaF) | ...
     (mask.deltaF  ~= deltaF) )
  error('deltaF mismatch');
end

% check that fRef is in range
fhigh = flow + (numFreqs-1)*deltaF;
if not(useHf)
  if ( (fRef < flow) | (fRef > fhigh) )
    error('fRef not in range');
  end
end
                                                                                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% frequency values
f = flow + transpose(deltaF*[0:1:numFreqs-1]);

% calculate normalization integral
[integral, integrand] = ...
  normalizationIntegral(gamma, fRef, alpha, calP1, calP2, mask);
% (Factor of 1/2 in heterodyned case taken care of in normalizationIntegral.m)

% calculate constant multiplicative and window factors 
H100 = HubbleConstant;
if useHf
  constfac = 1;
else
  constfac = (10 * pi^2)/(3 * H100^2);
end;

if window1 == -1 & window2 == -1
  w1w2bar = pp.w1w2bar;
  w1w2squaredbar = pp.w1w2squaredbar;
  w1w2ovlsquaredbar = pp.w1w2ovlsquaredbar;
else
  [w1w2bar, w1w2squaredbar, w1w2ovlsquaredbar]=windowFactors(window1, window2);
end

% calculate normalization and theoretical variance
norm = constfac * 2 * (1/w1w2bar) * (1/integral);
theorVar = T * (constfac^2) * (w1w2squaredbar/w1w2bar^2) * (1/integral);

% calculate sensitivity integrand
sensIntData = (1/T)*(1/constfac^2)*((w1w2bar^2)/w1w2squaredbar)*integrand.data;

% Define theorVar to be be the variance of the real part, in the case
% of heterodyned data.  Similarly for sensitivity integrand
if ( symmetry == 0 )
  theorVar = theorVar/2;
  sensIntData = 2*sensIntData;
end

% construct sensitivity integrand freq-series structure
sensInt = constructFreqSeries(sensIntData, integrand.flow, integrand.deltaF, ...
                              integrand.symmetry);

return

