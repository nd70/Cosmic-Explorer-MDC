function [integral, integrand] = ...
  normalizationIntegral(gamma, fRef, alpha, calP1, calP2, mask)
%
%  normalizationIntegral --- calculate the integral and integrand used
%  in the stochastic normalization
%
%  normalizationIntegral(gamma, fRef, alpha, calP1, calP2, mask)
%  calculates the integral and integrand
%
%  2 * int_flow^fhigh deltaF gamma^2(f) (f/fRef)^(2*alpha) /f^6 P1(f) P2(f) 
%
%  which is needed to evaluate the normalization constant and theoretical 
%  variance of the optimally-filtered CC statistic for a SB with power-law
%  frequency dependence.  The integrand is a frequency-series structure 
%  with associated metadata.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Extension for arbitrary signal power spectrum:
%
%  function [integral, integrand] = ...
%    normalizationIntegralRadiometer(gamma, Hf, [], calP1, calP2, mask)
%
%  i.e. instead of fref it takes a frequency series. This version is used
%  for the radiometer, but could also be usefull for a full-sky search.
%  It calculates the integral and integrand
%
%  2 * int_flow^fhigh deltaF abs(gamma)^2(f) Hf(f)^2 /  P1(f) P2(f) 
%
%
%  NOTE: mask zeroes out certain frequency bins
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org
%  Extendend by Stefan Ballmer, sballmer@ligo.mit.edu
%  
%  $Id: normalizationIntegral.m,v 1.12 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
useHf=true;
try
  fRef.flow;
catch
  useHf = false;
end;
if useHf
  % do error checking
  if ( (length(fRef.data) ~= numFreqs) )
    error('size mismatch');
  end
  if ( (fRef.flow ~= flow) )
    error('flow mismatch');
  end
  if ( (fRef.deltaF ~= deltaF) )
    error('deltaF mismatch');
  end  
end
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
f = flow + deltaF*transpose([0:1:numFreqs-1]);

% calculate integrand
if useHf
  numerator   = mask.data .* (conj(gamma.data).*gamma.data) .* (fRef.data.^2);
  denominator = calP1.data .* calP2.data;
else
  numerator   = mask.data .* (conj(gamma.data).*gamma.data) .* ((f/fRef).^(2*alpha));
  denominator = (f.^6) .* calP1.data .* calP2.data;
end
if ( flow == 0 )
  % ignore first terms in the arrays (set DC to zero)
  data = [ 0.0 ; 2 * numerator(2:end)./denominator(2:end) ];
else
  data = 2 * numerator./denominator;
end

% renormalize for heterodyned data
if ( symmetry == 0 )
  data = data / 2;
end

% construct frequency serires
integrand = constructFreqSeries(data, flow, deltaF, symmetry);

% calculate integral
integral = sum(integrand.data)*deltaF;

return
