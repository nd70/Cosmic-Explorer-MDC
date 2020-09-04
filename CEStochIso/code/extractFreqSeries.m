function [data, flow, deltaF, symmetry] = extractFreqSeries(x)
%
%  extractFreqSeries --- Extract frequency-series data from structure
%
%  extractFreqSeries(x) returns the data and metadata associated with 
%  a frequency-series data structure.
%
%  Outputs:
%
%  data      =  array of frequency-series values
%  flow      =  start frequency in Hz
%  deltaF    =  sampling rate in Hz
%  symmetry  =  the implied symmetry property of the frequency-series
%            =  1 if negative frequency values are the complex conjugates
%               of the positive frequency values (e.g., FFT of real data 
%               or 2-sided power spectra of real data)
%            =  0 if no implied symmetry (e.g., 1-sided power spectra or
%               FFT of complex data)
%
%  If the frequency-series symmetry metadata is missing, it is set to 0.
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org
%
%  $Id: extractFreqSeries.m,v 1.4 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = x.data;
flow = x.flow;
deltaF = x.deltaF;

try 
  symmetry = x.symmetry;
catch
  symmetry = 0;
end

return
