function [sensInts, flow, deltaF, gpsTimes] = readSensIntsFromFile(filename)
%
%  readSensIntsFromFile -- read sensitivity integrands from file
%
%  [sensInts, flow, deltaF, gpsTimes] = readSensIntsFromFile(filename)
%  returns a 2-d array of sensitivity integrands read in from a file 
%  (integrands corresponding to different times are in different rows).  
%  Also returned are the GPS times and initial frequency and frequency 
%  spacing corresponding to the sensitivity integrands.  The sens ints 
%  are sorted according to GPS time.
%
%  Routine written by Joseph D. Romano.
%  Contact Joseph.Romano@astro.cf.ac.uk
%
%  $Id: readSensIntsFromFile.m,v 1.1 2005-04-07 14:01:00 jromano Exp $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % default return values
  sensInts = [];
  flow = [];
  deltaF = [];
  gpsTimes = [];
  
  checkFileExists('Sensitivity integrand file', filename);
  
  [pathstr, name, ext] = fileparts(filename);
  if (ext == '.mat');

    % Check if there are no intervals to include, which can happen if eg. the science segment
    % was too short. If there are no segments we return
    load(filename, 'params');
    if (params.numIntervalsTotal < 1)
      return;
    end;

    load(filename, 'segmentStartTime', 'sensInt');
    gpsTimes = segmentStartTime;
    % Frequency series characteristics are the same for all spectra so we use
    % the first one to set them.
    flow  = sensInt(1, 1).flow;
    deltaF  = sensInt(1, 1).deltaF;
    numFreqs = length(sensInt(1,1).data);

    numGPSTimes = length(gpsTimes);
    sensInts = zeros(numGPSTimes, numFreqs);
    for k = 1:numGPSTimes
      sensInts(k, :) = sensInt(k, 1).data;
    end;
    clear segmentStartTime sensInt;
  else
    % Read in data from text file
    % For individual jobs files, there will usually be several spectra in a single file,
    % grouped into rows with the same GPS time.
    [gpsTimes, ccSigmas, freqs, sensInts] ...
        = textread(filename, '%f%f%f%f\n', -1, 'commentstyle', 'matlab');

    if isempty(gpsTimes)
      return;
    end
    
    % This pulls out the indices corresponding to each unique GPS time
    [gpsTimes, ind] = unique(gpsTimes);
    freqs = unique(freqs);
                                                                                
    numFreqs = length(freqs);
    flow = freqs(1);
    deltaF = freqs(2) - freqs(1);
    
    % create matrix with rows labeled by gps times and columns by frequencies
    sensInts = reshape(sensInts, numFreqs, length(gpsTimes));
    sensInts = transpose(sensInts);
  end;
    
  % sort the arrays on gps start times (if not already sorted)
  [gpsTimes, ind] = sort(gpsTimes);
  sensInts = sensInts(ind, :);

return;
