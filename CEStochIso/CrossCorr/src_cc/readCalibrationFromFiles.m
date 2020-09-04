function [t, f, R0, C0, alpha, gamma] = ...
  readCalibrationFromFiles(alphaBetaFile, calResponseFile, calCavGainFile)
%
%  readCalibrationFromFiles --- reads in calibration information from matfiles 
%
%  readCalibrationFromFiles(alphaBetaFile, calResponseFile, calCavGainFile)
%  reads in calibration information from matfiles, returning an array of
%  discrete times, discrete frequencies, values of the reference response
%  function (R0), values of the reference sensing function (C0), values of
%  the cavity factor alpha (C = alpha*C0), and values of the open loop factor 
%  gamma (H = gamma*[C0*R0 - 1]).
%
%  NOTE: alpha, gamma are not treated as time-series structures, since the 
%  discrete time-values may not be evenly spaced; similarly for the 
%  reference spectra R0, C0.
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: readCalibrationFromFiles.m,v 1.6 2005-05-05 14:47:50 sballmer Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(alphaBetaFile);
load(calResponseFile);
load(calCavGainFile);

% set NaN's occuring at f=0 in Response function to zero
if f(1)==0
  R0(1)=0;
end

return
