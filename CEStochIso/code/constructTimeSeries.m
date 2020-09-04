function x = constructTimeSeries(data, tlow, deltaT, fbase, phase)
%
%  constructTimeSeries -- Construct time-series structure
%
%  constructTimeSeries(data, tlow, deltaT, fbase, phase) returns a data
%  structure containing a time series and its associated metadata.
%
%  data   =  array of time-series values
%  tlow   =  start time in seconds (usually GPS time)
%  deltaT =  sampling time in seconds (i.e., 1/the sampling frequency)
%  fbase  =  the base frequency (in Hz) used in heterodyning the data
%  phase  =  initial phase (in radians) of the complex exponential
%            used in heterodyning
%
%  The last two arguments are optional, and if one or both is missing,
%  the corresponding field(s) are set to NaN.
% 
%  The convention used for heterodyning is, in the analog time domain,
%  heterodyneddata = originaldata * exp(-i*(2*pi*fbase*t+phase))
%  This is consistent with the convention defined in
%  http://www.lsc-group.phys.uwm.edu/daswg/docs/technical/T010095.pdf
%
%  Routine written by Joseph D. Romano and John T. Whelan.
%  Contact Joseph.Romano@astro.cf.ac.uk and/or john.whelan@ligo.org
%
%  $Id: constructTimeSeries.m,v 1.8 2006-04-14 21:05:42 whelan Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x.data = data;
x.tlow = tlow;
x.deltaT = deltaT;

try
  x.fbase = fbase;
catch
  x.fbase = NaN;
end;

try
  x.phase = phase;
catch
  x.phase = NaN;
end;

return
