function [ccStat, ccSpec] = calCSD(rbartilde1, rbartilde2, Q, ...
                                         response1, response2,tmax,NFFT)
%  
%  calCSD --- calculates the value and spectrum of the CSD
%
%  calCSD(rbartilde1, rbartilde2, Q, response1, response2)
%  calculates the value and spectrum of the standard optimally-filtered 
%  CSD given DFTs of the zero-padded uncalibrated data and the.  
%  response1,2 are response functions 
%  (strain/count) appropriate for the data segment being analysed.
%
%  If tmax is specified ccStat is a time series containing the inverse
%  fourier transform of ccSpec. It is the CSD statistic
%  for a set of time shifts smaller that tmax, spaced
%  by deltaT=1/(NFFT * deltaF). NFFT is optional.
%
%  Routine written by Vuk Mandic and Stefan Ballmer
% 
%  $Id: calCSD.m,v 1.1 2007-03-31 00:10:36 vmandic Exp $
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

if ( rbartilde1.symmetry == 1 & rbartilde2.symmetry == 1 )
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

% check that frequencies actually overlap
if ( rbartilde1.flow > ...
     rbartilde2.flow + rbartilde2.deltaF * (length(rbartilde2.data)-1) ...
     | rbartilde2.flow > ...
     rbartilde1.flow + rbartilde1.deltaF * (length(rbartilde1.data)-1) )
  error('frequency ranges do not overlap');
end;

% check that the rbartildes have the same deltaF
if ( rbartilde1.deltaF ~= rbartilde2.deltaF )
  error('input data deltaF mismatch');
end;

% check that reponse functions and Q have the same flow
if ( (response1.flow ~= flow ) | (response2.flow ~= flow ) )
  error('response function, optimal filter flow mismatch');
end;

% check that reponse functions and Q have the same deltaF
if ( (response1.deltaF ~= deltaF ) | (response2.deltaF ~= deltaF ) )
  error('response function, optimal filter deltaF mismatch');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% trim rbartildes so they refer to the same frequencies

offset = (rbartilde2.flow-rbartilde1.flow)*(1/rbartilde1.deltaF);

% check that offset in flows is an integer number of freq bins
if ( offset ~= floor(offset) )
  error('input data flow mismatch');
end;
if offset > 0
  rbartilde1.flow = rbartilde2.flow;
  rbartilde1.data = rbartilde1.data((1+offset):end);
elseif offset < 0
  rbartilde2.flow = rbartilde1.flow;
  rbartilde2.data = rbartilde2.data((1-offset):end);
end;

if length(rbartilde1.data) > length(rbartilde2.data)
  rbartilde1.data = rbartilde1.data(1:length(rbartilde2.data));
elseif length(rbartilde2.data) > length(rbartilde1.data)
  rbartilde2.data = rbartilde2.data(1:length(rbartilde1.data));
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% form product of the input data
data = conj(rbartilde1.data) .* rbartilde2.data;
rr   = constructFreqSeries(data, rbartilde1.flow, rbartilde1.deltaF);

% coarse grain the product to agree with optimal filter
% NOTE: index1, index2, frac1, frac2 are not needed for our analysis
[rrCG, index1, index2, frac1, frac2]=coarseGrain(rr, flow, deltaF, numFreqs);

% calibrate the coarse-grained product of the input data
data = conj(response1.data) .* rrCG.data .* response2.data;
ccSpec = constructFreqSeries(data, flow, deltaF);

% calculate the integral of CSD
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
