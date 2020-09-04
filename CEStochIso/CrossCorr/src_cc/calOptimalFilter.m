function [Q, theorVar, sensInt] = ... 
  calOptimalFilter(T, gamma, fRef, alphaExp, calPSD1, calPSD2, ...
                   window1, window2, mask, pp)
% 
%  calOptimalFilter --- calculate the optimal filter for the CC statistic
%
%  calOptimalFilter(T, gamma, fRef, alphaExp, calPSD1, calPSD2, window1, 
%  window2, mask, pp) calculates the (fully-calibrated) optimal filter 
%  function Q for the standard CC statistic, normalized so that the 
%  expected value of the statistic is OmegaRef * T.  Also returns the 
%  theoretical variance of the CC statistic, and the integrand of 
%  1/(theoretical variance)--the so-called sensitivity integrand.
%  The optimal filter and sensitivity integrand are frequency-series 
%  structures with corresponding metadata.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Extension for arbitrary signal power spectrum:
%
%  function [Q, theorVar, sensInt] = ... 
%    calOptimalFilter(T, gamma, Hf, [], calPSD1, calPSD2, ...
%                     window1, window2, mask, pp)
%
%  Calculates the (fully-calibrated) optimal filter function Q for
%  the standard CC statistic, normalized so that the expected value
%  of the statistic is eta * T for a point source with strength
%  eta*H(f) at the same location as specified in gamma.
%
%  Note:if window1 = window2 = -1, one needs to specify the input parameter
%  pp, which is a structure containing the already calculated window 
%  parameters.
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org 
%  Extendend by Stefan Ballmer, sballmer@ligo.mit.edu
% 
%  $Id: calOptimalFilter.m,v 1.16 2007-07-02 10:50:48 elr Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
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
                                                                                
% check that frequency series have the correct length
useHf = true;
try
  fRef.flow;
catch
  useHf = false;
end;

if ( (length(gamma.data) ~= numFreqs) | ...
     (length(calPSD1.data) ~= numFreqs) | ...
     (length(calPSD2.data) ~= numFreqs) | ...
     (length(mask.data) ~= numFreqs) )
  error('size mismatch');
end
                                                                                
% check that frequency series have the correct flow
if ( (gamma.flow ~= flow) | ...
     (calPSD1.flow ~= flow) | ...
     (calPSD2.flow ~= flow) | ...
     (mask.flow ~= flow) )
  error('flow mismatch');
end
                                                                                
% check that frequency series have the same deltaF
if ( (gamma.deltaF ~= deltaF) | ...
     (calPSD1.deltaF ~= deltaF) | ...
     (calPSD2.deltaF ~= deltaF) | ...
     (mask.deltaF ~= deltaF) )
  error('deltaF mismatch');
end
                                                                                % check that fRef is in range
fhigh = flow + (numFreqs-1)*deltaF;
if not(useHf)
  if ( (fRef < flow) | (fRef > fhigh) )
    error('fRef not in range');
  end 
end 

%check if windows are defined
if ~(all(window1 == -1) && all(window2 == -1))
  pp = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% frequency values
f = flow + deltaF*transpose([0:1:numFreqs-1]);

% calculate normalization and theoretical variance
[norm, theorVar, sensInt] = normalization(T, gamma, ...
                                          fRef, alphaExp, ...
                                          calPSD1, calPSD2, ...
                                          window1, window2, mask, pp);

% calculate optimal filter

if useHf
  numerator = norm * mask.data .* gamma.data .* fRef.data;
  denominator = calPSD1.data .* calPSD2.data;
else
  numerator = norm * mask.data .* gamma.data .* ( (f/fRef).^alphaExp );
  denominator = (f.^3) .* calPSD1.data .* calPSD2.data;
end

if (flow == 0)
  % ignore first terms in the arrays (set DC to zero)
  data = [ 0.0 ; numerator(2:end)./denominator(2:end) ]; 
else
  data = numerator./denominator;
end

% fill structure
Q = constructFreqSeries(data, flow, deltaF, symmetry);

return

