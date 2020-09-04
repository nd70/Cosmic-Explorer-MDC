function [ccStat, ccSpec] = processCSD(CSD, Q,tmax,NFFT)
%  
%  processCSD is samce as calCrossCorr, except that it takes the CSD  as input
%  calCrossCorr --- calculates the value and spectrum of the CC statistic
%
%  calCrossCorr(rbartilde1, rbartilde2, Q, response1, response2)
%  calculates the value and spectrum of the standard optimally-filtered 
%  CC statistic given DFTs of the zero-padded uncalibrated data and the 
%  calibrated optimal filter.  response1,2 are response functions 
%  (strain/count) appropriate for the data segment being analysed.
%
%  If tmax is specified ccStat is a time series containing the inverse
%  fourier transform of ccSpec. It is the CC statistic
%  for a set of time shifts smaller that tmax, spaced
%  by deltaT=1/(NFFT * deltaF). NFFT is optional.
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org
%  Extendend by Stefan Ballmer, sballmer@ligo.mit.edu
% 
%  $Id: processCSD.m,v 1.1 2007-04-04 19:07:40 vmandic Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try 
  tmax; 
catch 
  tmax = []; 
end;
try
  NFFT; 
catch 
  NFFT = []; 
end;
% extract frequency series metadata 
[data, flow, deltaF, symmetry] = extractFreqSeries(Q);
numFreqs = length(data);

heterodyned = (symmetry == 0);

if ( CSD.symmetry == 1)
  if heterodyned
    error('Trying to do heterodyned analysis on non-heterodyned data');
  end;
else
  if (~ heterodyned)
    error('Trying to do non-heterodyned analysis on heterodyned data');
  end;
end;

% error checking %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check that flow is >= 0
if ( flow < 0 )
  error('flow < 0');
end;

if ( flow == 0 & heterodyned )
  error('attempt to include DC for heterodyned data!')
end;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% form product of the input data
%data = conj(rbartilde1.data) .* rbartilde2.data;
%rr   = constructFreqSeries(data, rbartilde1.flow, rbartilde1.deltaF);

% coarse grain the product to agree with optimal filter
% NOTE: index1, index2, frac1, frac2 are not needed for our analysis
%[rrCG, index1, index2, frac1, frac2]=coarseGrain(rr, flow, deltaF, numFreqs);

% calibrate the coarse-grained product of the input data
%data = conj(response1.data) .* rrCG.data .* response2.data;
%ssCG = constructFreqSeries(data, flow, deltaF);

% calculate the integrand of the CC statistic 
data   = CSD.data.*Q.data;
ccSpec = constructFreqSeries(data, flow, deltaF);

% calculate the value of the CC statistic
if length(tmax)==0
  if flow == 0
    % dc component is real
    ccStat = ccSpec.data(1) + 2 * real(sum(ccSpec.data(2:end)));
  elseif heterodyned
    ccStat = sum(ccSpec.data);
  else
    ccStat = 2 * real(sum(ccSpec.data));
  end;
  ccStat = ccStat*deltaF;
else
  ccSpec.symmetry=1;
  ccStat=invFFT(ccSpec,tmax,NFFT);
end;

return
