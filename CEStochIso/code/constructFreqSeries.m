function x = constructFreqSeries(data, flow, deltaF, symmetry)
%
%  constructFreqSeries --- Construct frequency-series structure
%
%  constructFreqSeries(data, flow, deltaF, symmetry) returns a data
%  structure containing a frequency series and its associated metadata.
%
%  data     =  array of frequency-series values
%  flow     =  start frequency in Hz (e.g., for the FFT of heterodyned 
%              data, flow is the base frequency minus Nyquist)
%  deltaF   =  sampling rate in Hz
%  symmetry =  1 if negative frequency values are the complex conjugates
%              of the positive frequency values (e.g., FFT of real data or
%              2-sided power spectra of real data)  
%           =  0 if no implied symmetry (e.g., 1-sided power spectra or 
%              FFT of complex data)
%
%  The last argument is optional; if missing it is set to 0
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org
%
%  $Id: constructFreqSeries.m,v 1.4 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x.data = data;
x.flow = flow;
x.deltaF = deltaF;

try 
  x.symmetry = symmetry;
catch
  x.symmetry = 0;
end

return
